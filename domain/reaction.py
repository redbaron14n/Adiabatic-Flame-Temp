# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Reaction Class File
# ###################

from chempy import balance_stoichiometry
from domain.compound import Compound
from domain.compounds import compounds
from config import products_from_reactants
import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq

class Reaction:

    def __init__(
        self,
        reactants: set[str],
        temperatures: dict[str, float],
        dissociation: bool = False,
    ):  # Potentially arguments for reaction complexity
        """
        Initializes a Reaction object given a set of reactant Compounds.

        @param reactants : set[str] - Set of Compound.id strings representing the reactants of the reaction.
        @param temperatures : dict[str, float] - Dictionary mapping each reactant Compound.id to its entry temperature (K).
        @param dissociation : bool - Flag indicating whether to consider dissociation in the reaction (default is False).

        @attrib reactants : set[str] - Set of Compound objects representing the reactants of the reaction.
        @attrib inert_reactants : set[str | None] - Set of Compound objects representing inert reactants, if any
        @attrib products : set[str] - Set of Compound objects representing the products of the reaction.
        @attrib stoichiometry : tuple[dict[str, int], dict[str, int]] - Tuple containing two dictionaries representing the stoichiometric coefficients of reactants and products.
        @attrib delta_Hf : float - Total formation enthalpy change (kJ) for the reaction.
        """

        self._set_reactants(reactants)
        self._set_inert_reactants(dissociation)
        self._set_products(dissociation)
        self._set_stoichiometry()
        self._set_temperatures(temperatures)
        self._set_temperature_bounds()


    def _set_reactants(self, reactants: set[str]):

        self.reactants = reactants


    def _set_inert_reactants(self, dissociation: bool):

        self.inert_reactants = products_from_reactants(self.reactants, dissociation)[1]


    def _set_products(self, dissociation: bool):

        self.products = products_from_reactants(self.reactants, dissociation)[0]


    def _set_stoichiometry(self):

        reactive_species = self.reactants - self.inert_reactants
        reactant_strs = {compounds[r].formula for r in reactive_species}
        product_strs = {compounds[p].formula for p in self.products}
        balanced_reactants, balanced_products = balance_stoichiometry(
            reactant_strs, product_strs
        )
        for inert in self.inert_reactants:
            balanced_reactants[compounds[inert].formula] = 0
        self.stoichiometry = (balanced_reactants, balanced_products)


    def _set_temperatures(self, temperatures: dict[str, float]):

        if len(temperatures) != len(self.reactants):
            raise ValueError("Number of temperatures provided does not match number of reactants.")
        self.temperatures = temperatures


    """
    Using initial concentrations of reactants, determines and returns proportion of reactants used up
    """
    def _find_extent_of_reaction(self, concentrations: dict[str, float]) -> float:

        reactive_species = self.reactants - self.inert_reactants
        weighted_conc = {c: concentrations[c] / self.stoichiometry[0][compounds[c].formula] for c in reactive_species}
        extent = min(weighted_conc.values())
        return extent


    """
    Using intitial concentrations and extent of reaction, calculates final amounts of every species; reactants and products
    """
    def _compute_final_species_amounts(self, concentrations: dict[str, float], extent: float,) -> dict[str, float]:

        final_amounts = {}
        for reactant in self.reactants:
            initial_amount = concentrations[reactant]
            if reactant in self.inert_reactants:
                final_amounts[reactant] = initial_amount
            else:
                coef = self.stoichiometry[0][compounds[reactant].formula]
                consumed_amount = extent * coef
                final_amounts[reactant] = initial_amount - consumed_amount
        for product in self.products:
            final_amounts[product] = extent * self.stoichiometry[1][compounds[product].formula]
        return final_amounts


    def _calc_Hf(self, initial_amounts: dict[str, float], final_amounts: dict[str, float]) -> float:

        delta_Hf = 0.0
        for product in self.products:
            delta_Hf += final_amounts[product] * compounds[product].stdHf
        for reactant in self.reactants:
            delta_Hf -= initial_amounts[reactant] * compounds[reactant].stdHf
        return delta_Hf


    """
    Calculates and returns sensible heat of reactants at their introduction temperatures
    """
    def _calc_SH_reactants(self, initial_amounts: dict[str, float]) -> float:

        total_SH = 0.0
        for reactant in self.reactants:
            temp = self.temperatures[reactant]
            total_SH += initial_amounts[reactant] * compounds[reactant].SH(temp)
        return total_SH


    """
    Calculates sensible heat of products and leftover (unreacted) reactants at a given temperature.
    """
    def _calc_SH_products(self, final_amounts: dict[str, float], temperature: float) -> float:

        total_SH = 0.0
        for product in self.products:
            total_SH += final_amounts[product] * compounds[product].SH(temperature)
        for reactant in self.reactants:
            total_SH += final_amounts[reactant] * compounds[reactant].SH(temperature) # Adds sensible heats of any unreacted reactants, treating them as faux-products.
        return total_SH


    """
    Helper function for calc_flame_temp
    """
    def _energy_balance(self, temperature: float, initial_amounts: dict[str, float], final_amounts: dict[str, float]) -> float:

        residual = float(self._calc_SH_products(final_amounts, temperature)
                         - self._calc_SH_reactants(initial_amounts)
                         + self._calc_Hf(initial_amounts, final_amounts))
        return residual


    """
    Finds minimum and maximum shared temperatures from reactant data. Prevents extrapolation
    """
    def _set_temperature_bounds(self):

        min_temp = 0.0
        max_temp = np.inf
        for component in self.reactants.union(self.products):
            temperatures = compounds[component].get_temperatures()
            min_temp = max(min_temp, np.min(temperatures))
            max_temp = min(max_temp, np.max(temperatures))
        self.min_temp = min_temp
        self.max_temp = max_temp


    def _validate_concentrations(self, conc_dict: dict[str, float]) -> None:

        total = sum(conc_dict.values())
        if not np.isclose(total, 1.0):
            raise ValueError(f"Concentrations do not sum to 1.0. Sum: {total:.5f}")


    """
    Uses initial concentrations of reactants to find at what temperature the sensible heat of the products is equal to the sensible heat of reactants and heat of formation of reaction
    """
    def calc_flame_temp(self, concentrations: dict[str, float]) -> float:

        self._validate_concentrations(concentrations)
        extent = self._find_extent_of_reaction(concentrations)
        final_amounts = self._compute_final_species_amounts(concentrations, extent)
        if (self._energy_balance(self.min_temp, concentrations, final_amounts) * self._energy_balance(self.max_temp, concentrations, final_amounts) > 0):  # No root (flame temp) in bounds
            flame_temp = np.nan
        else:
            result = brentq(self._energy_balance, self.min_temp, self.max_temp, args=(concentrations, final_amounts))
            flame_temp = result[0] if isinstance(result, tuple) else result
        return flame_temp


    """
    Turns ratios into proportions
    """
    def _normalize(self, d: dict[str, float]) -> dict[str, float]:

        total = sum(d.values())
        return {k: v/total for k, v in d.items()}
    

    """
    Given a controlled reactant and its current concentration, scales the remaining reactants to appropriate concentrations with fixed relative ratios
    """
    def _scale_dependents(self, variable: str, x: float, ratios: dict[str, float]) -> dict[str, float]:

        """
        @param (str) variable : The compound id corresponding to the controlled reactant
        @param (float) x : The current concentration of the controlled reactant
        @param (dict[str, float]) ratios : The dictionary mapping the ratios of all of the reactants, including the variable
        """

        dependents = {k: v for k, v in ratios.items() if k != variable}
        dependents = self._normalize(dependents)
        leftover = 1.0 - x
        return {k: leftover * prop for k, prop in dependents.items()}
    

    """
    Generates initial concentrations of all reactants for every concentration of controlled reactant
    """
    def _generate_concentrations(self, variable: str, base_concs: dict[str, float], resolution: int = 100) -> list[dict[str, float]]:

        base_ratios = self._normalize(base_concs)
        delta_x = 1.0 / (resolution + 1)
        conc_list = []
        x_val = delta_x
        while x_val < 1.0:
            conc_dict = {variable: x_val}
            scaled_dependents = self._scale_dependents(variable, x_val, base_ratios)
            conc_dict.update(scaled_dependents)
            conc_list.append(conc_dict)
            x_val += delta_x
        return conc_list


    """
    Calculates the flame temperature data points as a function of the variable compound's concentration. Returns as an array
    """
    def calc_flame_table(self, variable_compound: str, base_concentrations: dict[str, float | int], resolution: int = 100) -> NDArray[np.float64]:

        concentration_dicts = self._generate_concentrations(
            variable_compound, base_concentrations, resolution
        )
        x_values = []
        flame_temps = []
        for conc_dict in concentration_dicts:
            x_values.append(conc_dict[variable_compound])
            flame_temp = self.calc_flame_temp(conc_dict)
            flame_temps.append(flame_temp)
        x_values, flame_temps = np.array(x_values), np.array(flame_temps)
        flame_table = np.stack((x_values, flame_temps))
        return flame_table
