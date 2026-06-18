# ###################
# Ian Janes
# Prof. Don Lipkin
# Adiabatic Flame Temperature Calculator
# ###################

from chempy import balance_stoichiometry
from config import determine_basic_products, INERTS
from domain.compounds import compounds, compounds_by_formula
from domain.dissociation import Dissociation
from math import log10
from numpy.typing import NDArray
from scipy.optimize import fsolve, least_squares, OptimizeResult
import numpy as np

PENALTY_FACTOR: float = 1000
MIN_LOG: float = -50
MAX_LOG: float = 5
MAX_TEMP: float = 6000.0
INIT_TEMP_GUESS: float = 3000.0

class CombustionReaction:

    def __init__(
            self,
            fuels: dict[str, float],
            oxidants: dict[str, float],
            temps: dict[str, float],
            pressure: float = 1e5,
            conc_res: int = 100
        ):

        """
        Initializes a new instance of the CombustionReaction class.

        :param dict[str, float] fuels: A dictionary mapping fuel compounds to their respective amounts.
        :param dict[str, float] oxidants: A dictionary mapping oxidant compounds to their respective amounts.
        :param dict[str, float] temps: A dictionary mapping compounds to their respective temperatures [K].
        :param float pressure: The pressure of the combustion reaction [Pa].
        :param int conc_res: The concentration resolution for the calculation.
        """
        
        self._set_fuels(fuels)
        self._set_oxidants(oxidants)
        self._set_present_atoms()
        self.pressure = pressure
        self.temperatures = temps
        self.concentration_resolution = conc_res
        self._set_conc_list()
        self._set_dependents()
        self._set_independents()
        self._set_item_indices()
        self._set_residual_indices()
        self._set_dissociations()
        self._set_stoich()
        self._set_bounds()


    ########################################
    # Getters and Setters
    ########################################


    def _set_fuels(self, fuels: dict[str, float]):

        for fuel, amount in fuels.items():
            if fuel not in compounds:
                raise ValueError("Invalid fuel compound.")
            if amount <= 0:
                raise ValueError("Fuel amounts must be positive values.")
        self._fuels = self._normalize_ratio(fuels)


    def _set_oxidants(self, oxidants: dict[str, float]):

        for oxidant, amount in oxidants.items():
            if oxidant not in compounds:
                raise ValueError("Invalid oxidant compound.")
            if amount <= 0:
                raise ValueError("Oxidant amounts must be positive values.")
        self._oxidants = self._normalize_ratio(oxidants)


    def _set_present_atoms(self):

        atoms: set[int] = set()
        for fuel in self._fuels:
            atoms.update(compounds[fuel].atomic_composition().keys())
        for oxidant in self._oxidants:
            atoms.update(compounds[oxidant].atomic_composition().keys())
        self._present_atoms: list[int] = sorted(atoms)

    
    @property
    def pressure(self):

        return self._pressure
    

    @pressure.setter
    def pressure(self, value: float):

        if value <= 0:
            raise ValueError("Pressure must be a positive value.")
        self._pressure = value


    @property
    def temperatures(self):

        return self._temperatures
    

    @temperatures.setter
    def temperatures(self, value: dict[str, float]):

        for temp in value.values():
            if temp < 0:
                raise ValueError("Temperature must be a non-negative value.")
        self._temperatures = value


    @property
    def concentration_resolution(self):

        return self._conc_res
    

    @concentration_resolution.setter
    def concentration_resolution(self, value: int):

        if value < 2:
            raise ValueError("Concentration resolution must be at least 2.")
        self._conc_res = value


    def _set_conc_list(self):

        denom = self._conc_res + 1 # Allows for desired resolution without 0 and 1 ratios (no flame)
        conc_list: list[dict[str, float]] = []
        for i in range(1, denom):
            fuel_ratio = i / denom
            conc_list.append({fuel: fuel_mix*fuel_ratio for fuel, fuel_mix in self._fuels.items()})
            conc_list[i-1].update({oxidant: oxi_mix*(1-fuel_ratio) for oxidant, oxi_mix in self._oxidants.items()})
        self._conc_list = conc_list


    def _set_dependents(self):

        if "Methane" in self._fuels:
            self._dependents = ["Methane", "Oxygen", "Water"]
            self._inv_dep_matr = np.linalg.inv(np.array([[4, 0, 2], [1, 0, 0], [0, 2, 1]]))
        elif "Hydrogen" in self._fuels:
            self._dependents = ["Hydrogen", "Oxygen"]
            self._inv_dep_matr = np.linalg.inv(np.array([[2, 0], [0, 2]]))
        else:
            raise ValueError("Invalid fuel compound.")
        

    def _set_independents(self):

        indeps: set[str] = set()
        actives = set(self._fuels.keys()).union(set(self._oxidants.keys()))
        self._reactants = actives
        basic_products = determine_basic_products(actives)
        self._basic_products = basic_products
        relevant_compounds = basic_products.union(actives)
        for compound in relevant_compounds:
            if compound not in self._dependents:
                indeps.add(compound)
            for dissociate in compounds[compound].dissociates:
                if dissociate not in self._dependents:
                    indeps.add(dissociate)
        indeps_list = sorted(indeps)
        self._independents = indeps_list
        self._set_indep_matr()


    def _set_indep_matr(self):

        shape = (len(self._present_atoms), len(self._independents))
        indep_matr: NDArray[np.float64] = np.zeros(shape)
        for i, indep in enumerate(self._independents):
            compound = compounds[indep]
            for j, atom in enumerate(self._present_atoms):
                indep_matr[j, i] = compound.atomic_composition().get(atom, 0)
        self._indep_matr = indep_matr


    def _set_item_indices(self):

        """
        Allows for consistent indexing of the guess array for the solver, if indexing is ever uncertain. More of a comfort thing.
        """

        self._item_indices = {item: idx for idx, item in enumerate(self._independents)}
        self._item_indices["T"] = len(self._independents) # Temperature is the last item in the guess array.
        offset = len(self._independents) + 1
        self._item_indices.update({dep: idx + offset for idx, dep in enumerate(self._dependents)})


    def _set_residual_indices(self):

        all_species = self._independents + self._dependents
        products: list[str] = []
        for compound in all_species:
            if compounds[compound].composition: # If not elemental
                products.append(compound)
        self._residual_indices = {product: idx for idx, product in enumerate(products)}
        self._residual_indices["T"] = len(self._residual_indices)


    def _set_dissociations(self):

        self._dissociations: dict[str, Dissociation] = {}
        for compound in self._residual_indices.keys():
            if compound == "T":
                continue
            components = set(compounds[compound].composition.keys())
            self._dissociations[compound] = Dissociation(compound, components)


    def _set_stoich(self):

        reactants = self._reactants
        products = self._basic_products
        inerts = reactants & INERTS
        reactants -= inerts
        reactant_frmls = {compounds[r].formula for r in reactants}
        product_frmls = {compounds[p].formula for p in products}
        dirty_react, dirty_prod = balance_stoichiometry(reactant_frmls, product_frmls)
        self._stoich: dict[str, float] = {}
        for r_frml, r_coeff in dirty_react.items():
            self._stoich[compounds_by_formula[r_frml].id] = int(r_coeff)
        for p_frml, p_coeff in dirty_prod.items():
            self._stoich[compounds_by_formula[p_frml].id] = int(p_coeff)
        for inert in inerts:
            self._stoich[inert] = 0


    def _set_bounds(self):

        guess_length = len(self._independents) + 1
        lower_log_bounds = np.full(guess_length - 1, MIN_LOG, dtype=np.float64)
        upper_log_bounds = np.full(guess_length - 1, MAX_LOG, dtype=np.float64)
        lower_temp_bound = np.array([0.0], dtype=np.float64)
        upper_temp_bound = np.array([MAX_TEMP], dtype=np.float64)
        lower_bounds = np.concatenate((lower_log_bounds, lower_temp_bound))
        upper_bounds = np.concatenate((upper_log_bounds, upper_temp_bound))
        self._bounds = (lower_bounds, upper_bounds)


    ########################################
    # Private Methods
    ########################################


    def _normalize_ratio(self, ratio: dict[str, float]) -> dict[str, float]:

        total = sum(ratio.values())
        return {compound: amount / total for compound, amount in ratio.items()}
    

    def _calc_init_atoms(self, conc: dict[str, float]) -> NDArray[np.float64]:

        """
        For given fuel and oxidant concentrations, calculates the initial number of each type of atom in the mixture.
        """

        init_atoms: NDArray[np.float64] = np.zeros(len(self._dependents))
        for compound_id, amount in conc.items():
            compound = compounds[compound_id]
            for i, atom in enumerate(self._present_atoms):
                init_atoms[i] += amount * compound.atomic_composition().get(atom, 0)
            # init_atoms[0] += amount * compound.atomic_composition().get(1, 0)
            # init_atoms[1] += amount * compound.atomic_composition().get(6, 0)
            # init_atoms[2] += amount * compound.atomic_composition().get(8, 0)
        return init_atoms
    

    def _calc_guess_vector(self, init_atoms: NDArray[np.float64], input_log_guess: NDArray[np.float64]) -> NDArray[np.float64]:

        """
        For a given guess of the independent species log quantities and temperature, calculates the full intended guess vector, including dependent species quantities.

        :param NDArray[np.float64] init_atoms: The initial number of each type of atom in the mixture.
        :param NDArray[np.float64] input_log_guess: The current guess of the independent species log quantities and temperature.
        :return NDArray[np.float64]: The calculated guess vector including dependent species quantities.
        """

        temp_index = self._residual_indices["T"]
        indep_qty_vector: NDArray[np.float64] = 10**np.delete(input_log_guess, temp_index)
        dep_qty_vector = self._inv_dep_matr @ (init_atoms - (self._indep_matr @ indep_qty_vector))
        full_real_guess = np.concatenate([indep_qty_vector, [input_log_guess[temp_index]], dep_qty_vector])
        return full_real_guess
    

    def _calc_extnt_of_react(self, conc_dict: dict[str, float]) -> float:

        reactives = self._reactants - INERTS
        weighted_conc = {c: conc_dict[c] / self._stoich[c] for c in reactives}
        return min(weighted_conc.values())
    

    def _calc_basic_final_amouns(self, conc_dict: dict[str, float], extent: float) -> dict[str, float]:

        """
        Calculates the final amounts of each compound after reaction if no dissociation were to occur. Used to find a reasonable initial guess.
        """

        final_amounts: dict[str, float] = {}
        inerts = self._reactants & INERTS
        reactants = self._reactants - inerts
        for r in reactants:
            coef = self._stoich[r]
            consumed = extent * coef
            final_amounts[r] = conc_dict[r] - consumed
        for i in inerts:
            final_amounts[i] = conc_dict[i]
        for p in self._basic_products:
            coef = self._stoich[p]
            final_amounts[p] = extent * coef
        return final_amounts
    

    def _calc_init_log_guess(self, conc_dict: dict[str, float]) -> NDArray[np.float64]:

        """
        Calculates the initial log-space guess for the equilibrium calculation. Trimmed to only independents.
        """

        extent = self._calc_extnt_of_react(conc_dict)
        basic_final_amounts = self._calc_basic_final_amouns(conc_dict, extent)
        init_log_guess: NDArray[np.float64] = np.full(len(self._independents) + 1, MIN_LOG, dtype=np.float64)
        for item in (basic_final_amounts.keys() & self._independents):
            if basic_final_amounts[item] == 0:
                init_log_guess[self._item_indices[item]] = MIN_LOG
            else:
                init_log_guess[self._item_indices[item]] = log10(basic_final_amounts[item])
        init_log_guess[self._item_indices["T"]] = INIT_TEMP_GUESS
        return init_log_guess
    

    def _convert_guess_to_log(self, full_real_guess: NDArray[np.float64]) -> NDArray[np.float64]:

        full_log_guess: NDArray[np.float64] = np.full(len(full_real_guess), MIN_LOG, dtype=np.float64)
        for item, i in self._item_indices.items():
            real_val = full_real_guess[i]
            if item == "T":
                full_log_guess[i] = real_val # The only value that is inputted in real-space.
                continue
            if real_val == 0:
                full_log_guess[i] = MIN_LOG
                continue
            full_log_guess[i] = log10(real_val)
        return full_log_guess
    

    def _calc_reactant_heat(self, conc_dict: dict[str, float]) -> float:

        reactant_heat: float = 0.0
        for r in self._reactants:
            qty = conc_dict[r] # Conc_dicts are in real-space
            react = compounds[r]
            entry_temp = self._temperatures[r]
            reactant_heat += qty*(react.SH(entry_temp) + react.stdHf)
        return reactant_heat
    

    def _calc_product_heat(self, full_log_guess: NDArray[np.float64]) -> float:

        temp = full_log_guess[self._item_indices["T"]]
        product_heat: float = 0.0
        for p, i in self._item_indices.items():
            if p == "T":
                continue
            qty = 10**full_log_guess[i]
            prod = compounds[p]
            product_heat += qty*(prod.SH(temp) + prod.stdHf)
        return product_heat
    

    def _calc_energy_residual(self, full_log_guess: NDArray[np.float64], conc_dict: dict[str, float]) -> float:

        reactant_heat = self._calc_reactant_heat(conc_dict)
        product_heat = self._calc_product_heat(full_log_guess)
        return product_heat - reactant_heat


    def _residual_function(self, init_log_guess: NDArray[np.float64], init_atoms: NDArray[np.float64], conc_dict: dict[str, float]) -> NDArray[np.float64]:

        full_real_guess = self._calc_guess_vector(init_atoms, init_log_guess)
        # print(f"Full real guess: {full_real_guess}")
        residual_penalty = 0.0
        if np.any(full_real_guess < -1e-10): # Set all residuals to a function of the OoB value method
            residual_penalty = PENALTY_FACTOR * (-np.min(full_real_guess))
        full_log_guess = self._convert_guess_to_log(np.maximum(full_real_guess, 10**MIN_LOG))
        print(full_log_guess)
        residuals = np.full(len(init_log_guess), 0.0, dtype=np.float64)
        for comp, diss_obj in self._dissociations.items():
            residuals[self._residual_indices[comp]] = diss_obj.equilibrium_residual(full_log_guess, self._item_indices)
        residuals[self._residual_indices["T"]] = self._calc_energy_residual(full_log_guess, conc_dict)
        return residuals + residual_penalty


    ########################################
    # Public Methods
    ########################################


    def equilibrate(self, conc_dict: dict[str, float]):

        init_log_guess = self._calc_init_log_guess(conc_dict)
        init_atoms = self._calc_init_atoms(conc_dict)
        equil = least_squares(self._residual_function, init_log_guess, args=(init_atoms, conc_dict), bounds=self._bounds, xtol=1e-12)
        return equil



    def calculate_temperatures(self):

        temps: list[float] = []
        for conc_dict in self._conc_list:
            equil = self.equilibrate(conc_dict)