# ###################
# Ian Janes
# Prof. Don Lipkin
# Adiabatic Flame Temperature Calculator
# ###################

from config import determine_basic_products
from domain.compounds import compounds
from domain.dissociation import Dissociation
from numpy.typing import NDArray
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


    def _set_dissociations(self):

        self._dissociations: dict[str, Dissociation] = {}
        for compound in self._residual_indices.keys():
            components = set(compounds[compound].composition.keys())
            self._dissociations[compound] = Dissociation(compound, components)


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
        :param NDArray[np.float64] input_guess: The current guess of the independent species quantities and temperature, in log form for the species quantities.
        :return NDArray[np.float64]: The calculated guess vector including dependent species quantities.
        :return bool: Whether the calculated dependent quantities are physically valid (non-negative).
        """

        indep_qty_vector = np.zeros(len(self._independents))
        for i, indep in enumerate(self._independents):
            indep_qty_vector[i] = 10**input_guess[i] # Convert log quantities to actual quantities.
        dep_qty_vector = self._inv_dep_matr @ (init_atoms - (self._indep_matr @ indep_qty_vector))
        temp_index = self._item_indices["T"]
        true_guess = np.concatenate([indep_qty_vector, [input_guess[temp_index]], dep_qty_vector])
        is_valid = bool(np.all(true_guess >= 0))
        return true_guess, is_valid
    

    def _calc_residuals(self, guess: NDArray[np.float64]) -> NDArray[np.float64]:

        residuals = np.full(len(self._independents), 1000, dtype=np.float64)
        for diss_reaction in self._dissociations.values():
            resid = diss_reaction.equilibrium_residual(guess, self._item_indices, self._pressure)
            indx = self._residual_indices[diss_reaction.molecule_id]
            residuals[indx] = resid
        return residuals
    

    def _residual_function(self, conc_idx: int, guess: NDArray[np.float64]) -> NDArray[np.float64]:

        """
        Process input arguments and calculates the residuals for a given concentration index and guess vector.
        """

        conc_dict = self._conc_list[conc_idx]
        print(conc_dict)
        init_atoms = self._calc_init_atoms(conc_dict)
        true_guess, is_valid = self._calc_guess_vector(init_atoms, guess)
        print(true_guess)
        residual = self._calc_residuals(true_guess)
        print(residual)
        # if not is_valid: # Blow up a non-physical solution
        #     min_mag = -np.min(true_guess) # Always positive
        #     residual *= 10**(NPS_EXP_FACTOR * min_mag)
        return residual