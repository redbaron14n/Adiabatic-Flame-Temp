# ###################
# Ian Janes
# Prof. Don Lipkin
# Adiabatic Flame Temperature Calculator
# ###################

from config import determine_basic_products
from domain.compounds import compounds
from domain.dissociation import Dissociation
from numpy.typing import NDArray
from scipy.optimize import fsolve, least_squares, OptimizeResult
import numpy as np

NPS_EXP_FACTOR: float = 3000 # The factor by which to multiply the exponent of the residual of the non-physical solution. Enforces mass balance. Arbitrary number

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
        self.pressure = pressure
        self.temperatures = temps
        self.concentration_resolution = conc_res
        self._set_conc_list()
        self._set_dependents()
        self._set_independents()
        self._set_item_indices()
        self._set_residual_indices()
        self._set_dissociations()
        # self._set_bounds()


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
            self._dependents = ["Hydrogen", "Oxygen", "Water"]
            self._inv_dep_matr = np.linalg.inv(np.array([[2, 0, 2], [0, 0, 0], [0, 2, 1]]))
        else:
            raise ValueError("Invalid fuel compound.")
        

    def _set_independents(self):

        indeps: set[str] = set()
        actives = set(self._fuels.keys()).union(set(self._oxidants.keys()))
        basic_products = determine_basic_products(actives)
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

        shape = (3, len(self._independents))
        indep_matr: NDArray[np.float64] = np.zeros(shape)
        for i, indep in enumerate(self._independents):
            compound = compounds[indep]
            indep_matr[0, i] = compound.atomic_composition().get(1, 0)
            indep_matr[1, i] = compound.atomic_composition().get(6, 0)
            indep_matr[2, i] = compound.atomic_composition().get(8, 0)
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


    # def _set_bounds(self):

    #     guess_length = len(self._independents)
    #     self._lower_bounds: NDArray[np.float64] = np.full(guess_length + 1, 0.0)
    #     self._upper_bounds: NDArray[np.float64] = np.concatenate((np.full(guess_length, 5.0), np.array([6000.0])))


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
            init_atoms[0] += amount * compound.atomic_composition().get(1, 0)
            init_atoms[1] += amount * compound.atomic_composition().get(6, 0)
            init_atoms[2] += amount * compound.atomic_composition().get(8, 0)
        return init_atoms
    

    def _calc_guess_vector(self, init_atoms: NDArray[np.float64], input_guess: NDArray[np.float64]) -> tuple[NDArray[np.float64], bool]:

        """
        For a given guess of the independent species log quantities and temperature, calculates the full intended guess vector, including dependent species quantities.

        :param NDArray[np.float64] init_atoms: The initial number of each type of atom in the mixture.
        :param NDArray[np.float64] input_guess: The current guess of the independent species log quantities and temperature.
        :return NDArray[np.float64]: The calculated guess vector including dependent species quantities.
        :return bool: Whether the calculated dependent quantities are physically valid (non-negative).
        """

        temp_index = self._residual_indices["T"]
        indep_qty_vector: NDArray[np.float64] = 10**np.delete(input_guess, temp_index)
        dep_qty_vector = self._inv_dep_matr @ (init_atoms - (self._indep_matr @ indep_qty_vector))
        full_guess = np.concatenate([indep_qty_vector, [input_guess[temp_index]], dep_qty_vector])
        is_valid = bool(np.all(full_guess >= 0))
        return full_guess, is_valid
    

    def _calc_equil_residual(self, compound: str, full_guess: NDArray[np.float64], is_valid: bool) -> float:

        if not is_valid:
            penalty = np.min(full_guess) * 1e4 # Need to go back to log form to make this actually potent
            return penalty
        residual = self._dissociations[compound].equilibrium_residual(full_guess, self._item_indices, self._pressure)
        return residual
    

    def _residual_function(self, guess, conc_dict):

        init_atoms = self._calc_init_atoms(conc_dict)
        full_guess, is_valid = self._calc_guess_vector(init_atoms, guess)
        residuals: NDArray[np.float64] = np.zeros(len(self._residual_indices))
        for compound in self._independents:
            residual = self._calc_equil_residual(compound, full_guess, is_valid)
            residuals[self._residual_indices[compound]] = residual
        residuals[self._residual_indices["T"]] = self._calc_heat_residual(full_guess, is_valid)
        return residuals
    

    ########################################
    # Public Methods
    ########################################


    def calculate_temperatures(self):

        temps: list[float] = []
        for conc_dict in self._conc_list:
            init_guess = self._calc_init_guess(conc_dict)
            equil = least_squares(self._residual_function, init_guess, args=(conc_dict,), bounds=self._bounds)