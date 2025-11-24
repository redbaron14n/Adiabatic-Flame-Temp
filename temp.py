# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Thermochemical Data File Editor
# ###################

from chempy import balance_stoichiometry
import numpy as np
from numpy.typing import NDArray
import pandas as pd
from scipy.interpolate import BSpline, make_interp_spline
from scipy.optimize import brentq

DATA_FILE = "thermochemical_data.csv"
TD = pd.read_csv(DATA_FILE)

def get_finite_list(list: NDArray) -> NDArray:

    finite_list = np.copy(list)
    finite_mask = np.isfinite(list)
    max_finite = np.max(list[finite_mask])
    finite_list[~finite_mask] = max_finite * 1e6
    return finite_list

def get_numerical(list: NDArray) -> tuple[NDArray, NDArray]:

    numerical_list = np.copy(list)
    for i in range(len(list)):
        if list[i] == "inf":
            numerical_list[i] = np.inf
    finite_list = get_finite_list(numerical_list)
    return numerical_list, finite_list

class Compound:

    ##### Initialization Methods #####

    def __init__(self, name: str, formula: str, id: str, ref_temp: float = 298.15):

        self.__set_name(name)
        self.__set_formula(formula)
        self.__set_id(id)
        self.__set_data(id)
        self.__set_ref_temp(ref_temp)

    def __set_name(self, name: str):

        self.name = name

    def __set_formula(self, formula: str):

        self.formula = formula

    def __set_id(self, id: str):

        self.__id = id

    def __set_data(self, id: str):

        self.data_table = TD[TD["Compound"] == self.__id]
        self.__set_temperatures(self.data_table["T"])
        self.__set_SH(self.data_table["SH"])
        self.__set_Hf(self.data_table["Hf"])
        self.__set_logKf(self.data_table["logKf"])

    def __set_temperatures(self, list: pd.Series):

        temp_array = list.to_numpy(dtype=np.float64)
        self.__temperatures = temp_array

    def __set_SH(self, list: pd.Series):

        SH_array = list.to_numpy()
        self.__SH_list = SH_array
        self.__SH_function = make_interp_spline(self.get_temperatures(), SH_array, k=1)

    def __set_Hf(self, list: pd.Series):

        Hf_array = list.to_numpy()
        self.__Hf_list = Hf_array
        self.__Hf_function = make_interp_spline(self.get_temperatures(), Hf_array, k=1)

    def __set_logKf(self, list: pd.Series):
        
        logKf_array = list.to_numpy()
        logKf_array, logKf_finite_array = get_numerical(logKf_array)
        self.__logKf_list = logKf_array
        self.__logKf_function = make_interp_spline(self.get_temperatures(), logKf_finite_array, k=1)

    def __set_ref_temp(self, ref_temp: float):

        self.ref_temp = ref_temp
        self.stdHf = self.Hf(ref_temp)

    ##### Attribute Methods #####

    def get_temperatures(self) -> NDArray[np.float64]:

        temp_array = self.__temperatures
        return temp_array
    
    def get_SH_list(self) -> NDArray[np.float64]:

        SH_array = self.__SH_list
        return SH_array
    
    def get_Hf_list(self) -> NDArray[np.float64]:

        Hf_array = self.__Hf_list
        return Hf_array
    
    def get_logKf_list(self) -> NDArray[np.float64]:

        logKf_list = self.__logKf_list
        return logKf_list
    
    ##### Functional Methods #####

    def SH(self, temperature: float) -> float:

        """
        Returns the sensible heat (kJ/mol) of the compound at a given temperature (K).
        """

        value = float(self.__SH_function(temperature))
        return value
    
    def Hf(self, temperature: float) -> float:

        """
        Returns the heat of formation (kJ/mol) of the compound at a given temperature (K).
        """

        value = float(self.__Hf_function(temperature))
        return value
    
    def logKf(self, temperature: float) -> float:

        """
        Returns the logKf of the compound at a given temperature (K).
        """

        value = self.__logKf_function(temperature)
        return value
    
class Reaction:

    def __init__(self, reactants: set[Compound], temperatures: dict[Compound, float], dissociation: bool = False): # Potentially arguments for reaction complexity

        """
        Initializes a Reaction object given a set of reactant Compounds.

        @param reactants : set[Compound] - Set of Compound objects representing the reactants of the reaction.
        @param temperatures : dict[Compound, float] - Dictionary mapping each reactant Compound to its entry temperature (K).
        @param dissociation : bool - Flag indicating whether to consider dissociation in the reaction (default is False).

        @attrib reactants : set[Compound] - Set of Compound objects representing the reactants of the reaction.
        @attrib products : set[Compound] - Set of Compound objects representing the products of the reaction.
        @attrib stoichiometry : tuple[dict[str, int], dict[str, int]] - Tuple containing two dictionaries representing the stoichiometric coefficients of reactants and products.
        @attrib delta_Hf : float - Total formation enthalpy change (kJ) for the reaction.
        """

        self.__set_reactants(reactants)
        self.__set_products(dissociation)
        self.__set_stoichiometry()
        self.__set_temperatures(temperatures)

    def __set_reactants(self, reactants: set[Compound]):

        self.reactants = reactants

    def __set_products(self, dissociation: bool):

        self.products = products_from_reactants(self.reactants, dissociation)

    def __set_stoichiometry(self):

        reactant_strs = {r.formula for r in self.reactants}
        product_strs = {p.formula for p in self.products}
        balanced_reactants, balanced_products = balance_stoichiometry(reactant_strs, product_strs)
        self.stoichiometry = (balanced_reactants, balanced_products)

    def __set_temperatures(self, temperatures: dict[Compound, float]):

        if len(temperatures) != len(self.reactants):
            raise ValueError("Number of temperatures provided does not match number of reactants.")
        self.temperatures = temperatures

    # def __calc_SH_of_reactants(self):

    #     """
    #     Calculates the total sensible heat (kJ) of the reactants at their respective entry temperatures.
    #     """

    #     total_SH = 0.0
    #     for reactant in self.reactants:
    #         coeff = self.stoichiometry[0][reactant.formula]
    #         temp = self.temperatures[reactant]
    #         total_SH += coeff * reactant.SH(temp)
    #     self.__total_SH_reactants = total_SH

    # def __calc_Hf(self):

    #     """
    #     Calculates the total formation enthalpy change (kJ) for the reaction based on the heats of formation of reactants and products at their respective temperatures.
    #     """

    #     delta_Hf = 0.0
    #     for reactant in self.reactants:
    #         coeff = self.stoichiometry[0][reactant.formula]
    #         delta_Hf -= coeff * reactant.stdHf
    #     for product in self.products:
    #         coeff = self.stoichiometry[1][product.formula]
    #         delta_Hf += coeff * product.stdHf
    #     self.delta_Hf = delta_Hf

    # def __create_SH_products_function(self):

    #     common_temps = set()
    #     for product in self.products:
    #         temps = [round(t, 6) for t in product.get_temperatures()]
    #         common_temps.update(temps)
    #     common_temps = sorted(common_temps)
    #     SH_functions = []
    #     stoich_products = self.stoichiometry[1]
    #     for product in self.products:
    #         stoich_coeff = stoich_products[product.formula]
    #         SH_vals = stoich_coeff * np.array([product.SH(t) for t in common_temps])
    #         function = make_interp_spline(common_temps, SH_vals, k=1)
    #         SH_functions.append(function)
    #     spline_coeffs = sum(f.c for f in SH_functions)
    #     self.SH_products_function = BSpline(SH_functions[0].t, spline_coeffs, k=1)

    def __validate_concentrations(self, concentrations: dict[Compound, float]):

        if concentrations.keys() != self.reactants:
            raise ValueError("Concentration keys do not match reactants.")
        elif not np.isclose(sum(concentrations.values()), 1.0):
            raise ValueError("Concentrations must sum to 1.")
        elif any(c <= 0 for c in concentrations.values()):
            raise ValueError("Concentrations must be greater than 0")
        elif any(c >= 1 for c in concentrations.values()):
            raise ValueError("Individual concentrations must be less than 1.")
    
    def __find_extent_of_reaction(self, concentrations: dict[Compound, float]) -> float:

        weighted_conc = {c: concentrations[c] / self.stoichiometry[0][c.formula] for c in self.reactants}
        extent = min(weighted_conc.values())
        return extent
    
    def __compute_final_species_amounts(self, concentrations: dict[Compound, float], extent: float) -> dict[Compound, float]:

        final_amounts = {}
        for reactant in self.reactants:
            initial_amount = concentrations[reactant]
            consumed_amount = extent * self.stoichiometry[0][reactant.formula]
            final_amounts[reactant] = initial_amount - consumed_amount
        for product in self.products:
            formed_amount = extent * self.stoichiometry[1][product.formula]
            final_amounts[product] = formed_amount
        return final_amounts
    
    def __calc_Hf(self, final_amounts: dict[Compound, float]) -> float:

        delta_Hf = 0.0
        for product in self.products:
            delta_Hf += final_amounts[product] * product.stdHf
        for reactant in self.reactants:
            delta_Hf -= final_amounts[reactant] * reactant.stdHf
        return delta_Hf
    
    def __calc_SH_reactants(self, final_amounts: dict[Compound, float]) -> float:

        total_SH = 0.0
        for reactant in self.reactants:
            temp = self.temperatures[reactant]
            total_SH += final_amounts[reactant] * reactant.SH(temp)
        return total_SH
    
    def __calc_SH_products(self, final_amounts: dict[Compound, float], temperature: float) -> float:

        total_SH = 0.0
        for product in self.products:
            total_SH += final_amounts[product] * product.SH(temperature)
        return total_SH
    
    def __energy_balance(self, temperature: float, final_amounts: dict[Compound, float]) -> float:

        residual = self.__calc_SH_products(final_amounts, temperature) - self.__calc_SH_reactants(final_amounts) + self.__calc_Hf(final_amounts)
        return residual
    
    def __get_temperature_bounds(self,) -> tuple[float, float]:

        min_temp = 0.0
        max_temp = np.inf
        for component in self.reactants.union(self.products):
            min_temp = max(min_temp, np.min(component.get_temperatures()))
            max_temp = min(max_temp, np.max(component.get_temperatures()))
        return min_temp, max_temp
    
    def calc_flame_temp(self, concentrations: dict[Compound, float]) -> float:

        self.__validate_concentrations(concentrations)
        extent = self.__find_extent_of_reaction(concentrations)
        final_amounts = self.__compute_final_species_amounts(concentrations, extent)
        T_min, T_max = self.__get_temperature_bounds()
        result = brentq(self.__energy_balance, T_min, T_max, args=(final_amounts,))
        flame_temp = result[0] if isinstance(result, tuple) else result
        return flame_temp

def products_from_reactants(reactants: set[Compound], dissociation: bool) -> set[Compound]: # Placeholder function for potential reaction product generation

    if (reactants == {Methane, Oxygen}) and not dissociation:
        products = {CarbonDioxide, Water}
    else:
        raise NotImplementedError("Reaction product generation not implemented for given reactants.")
    return products

CarbonDioxide = Compound("Carbon Dioxide", "CO2", "Carbon_Dioxide")
Methane = Compound("Methane", "CH4", "Methane")
Oxygen = Compound("Oxygen", "O2", "Oxygen")
Water = Compound("Water", "H2O", "Water")

test_reaction = Reaction({Methane, Oxygen}, {Methane: 300.0, Oxygen: 300.0}, dissociation=False)
flame_temp = test_reaction.calc_flame_temp({Methane: 0.2, Oxygen: 0.8})
print(flame_temp)