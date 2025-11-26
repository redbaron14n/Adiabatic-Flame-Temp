# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Reaction Class File
# ###################

from chempy import balance_stoichiometry
from compound import Compound
from config import products_from_reactants
import numpy as np
from numpy.typing import NDArray
import pandas as pd
from scipy.interpolate import make_interp_spline
from scipy.optimize import brentq
    
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

        self._set_reactants(reactants)
        self._set_products(dissociation)
        self._set_stoichiometry()
        self._set_temperatures(temperatures)
        self._set_temperature_bounds()

    def _set_reactants(self, reactants: set[Compound]):

        self.reactants = reactants

    def _set_products(self, dissociation: bool):

        self.products = products_from_reactants(self.reactants, dissociation)

    def _set_stoichiometry(self):

        reactant_strs = {r.formula for r in self.reactants}
        product_strs = {p.formula for p in self.products}
        balanced_reactants, balanced_products = balance_stoichiometry(reactant_strs, product_strs)
        self.stoichiometry = (balanced_reactants, balanced_products)

    def _set_temperatures(self, temperatures: dict[Compound, float]):

        if len(temperatures) != len(self.reactants):
            raise ValueError("Number of temperatures provided does not match number of reactants.")
        self.temperatures = temperatures

    def _validate_reactants(self, compounds: set[Compound]):

        reactant_set = set(self.reactants)
        if compounds != reactant_set:
            missing = reactant_set - compounds
            extra = compounds - reactant_set
            raise ValueError(f"Provided compounds do not match reactants. Missing: {missing}, Extra: {extra}")

    def _validate_concentrations(self, concentrations: dict[Compound, float]):

        self._validate_reactants(set(concentrations.keys()))
        if not np.isclose(sum(concentrations.values()), 1.0):
            raise ValueError("Concentrations must sum to 1.")
        elif any(c <= 0 for c in concentrations.values()):
            raise ValueError("Concentrations must be greater than 0")
        elif (len(concentrations) != 1) and any(c >= 1 for c in concentrations.values()):
            raise ValueError("Individual concentrations must be less than 1.")
    
    def _find_extent_of_reaction(self, concentrations: dict[Compound, float]) -> float:

        weighted_conc = {c: concentrations[c] / self.stoichiometry[0][c.formula] for c in self.reactants}
        extent = min(weighted_conc.values())
        return extent
    
    def _compute_final_species_amounts(self, concentrations: dict[Compound, float], extent: float) -> dict[Compound, float]:

        final_amounts = {}
        for reactant in self.reactants:
            initial_amount = concentrations[reactant]
            consumed_amount = extent * self.stoichiometry[0][reactant.formula]
            final_amounts[reactant] = initial_amount - consumed_amount
        for product in self.products:
            formed_amount = extent * self.stoichiometry[1][product.formula]
            final_amounts[product] = formed_amount
        return final_amounts
    
    def _calc_Hf(self, final_amounts: dict[Compound, float]) -> float:

        delta_Hf = 0.0
        for product in self.products:
            delta_Hf += final_amounts[product] * product.stdHf
        for reactant in self.reactants:
            delta_Hf -= final_amounts[reactant] * reactant.stdHf
        return delta_Hf
    
    def _calc_SH_reactants(self, final_amounts: dict[Compound, float]) -> float:

        total_SH = 0.0
        for reactant in self.reactants:
            temp = self.temperatures[reactant]
            total_SH += final_amounts[reactant] * reactant.SH(temp)
        return total_SH
    
    def _calc_SH_products(self, final_amounts: dict[Compound, float], temperature: float) -> float:

        total_SH = 0.0
        for product in self.products:
            total_SH += final_amounts[product] * product.SH(temperature)
        return total_SH
    
    def _energy_balance(self, temperature: float, final_amounts: dict[Compound, float]) -> float:

        residual = float(self._calc_SH_products(final_amounts, temperature) - self._calc_SH_reactants(final_amounts) + self._calc_Hf(final_amounts))
        return residual
    
    def _set_temperature_bounds(self):

        min_temp = 0.0
        max_temp = np.inf
        for component in self.reactants.union(self.products):
            min_temp = max(min_temp, np.min(component.get_temperatures()))
            max_temp = min(max_temp, np.max(component.get_temperatures()))
        self._min_temp = min_temp
        self._max_temp = max_temp
    
    def calc_flame_temp(self, concentrations: dict[Compound, float]) -> float:

        self._validate_concentrations(concentrations)
        extent = self._find_extent_of_reaction(concentrations)
        final_amounts = self._compute_final_species_amounts(concentrations, extent)
        if self._energy_balance(self._min_temp, final_amounts) * self._energy_balance(self._max_temp, final_amounts) > 0: # No root in bounds (flame temp higher than max)
            flame_temp = np.nan
        else:
            result = brentq(self._energy_balance, self._min_temp, self._max_temp, args=(final_amounts,))
            flame_temp = result[0] if isinstance(result, tuple) else result
        return flame_temp
    
    def _generate_concentrations(self, variable_compound: Compound, base_concentrations: dict[Compound, float], resolution: int = 100) -> list[dict[Compound, float]]:
        
        self._validate_concentrations(base_concentrations)
        delta_x = 1.0 / (resolution + 1) # Avoid 0 and 1 concentrations while maintaining resolution
        dependent_compounds = self.reactants - {variable_compound}
        total_dependent_conc = sum(base_concentrations[c] for c in dependent_compounds)
        concentration_list = []
        x_val = delta_x
        while x_val < 1.0:
            conc_dict = {variable_compound: x_val}
            for compound in dependent_compounds:
                base_conc = base_concentrations[compound]
                adjusted_conc = base_conc * (1.0 - x_val) / total_dependent_conc
                conc_dict[compound] = adjusted_conc
            concentration_list.append(conc_dict)
            x_val += delta_x
        return concentration_list
    
    def calc_flame_table(self, variable_compound: Compound, base_concentrations: dict[Compound, float], resolution: int = 100) -> NDArray[np.float64]:

        concentration_dicts = self._generate_concentrations(variable_compound, base_concentrations, resolution)
        x_values = []
        flame_temps = []
        for conc_dict in concentration_dicts:
            x_values.append(conc_dict[variable_compound])
            flame_temp = self.calc_flame_temp(conc_dict)
            flame_temps.append(flame_temp)
        x_values, flame_temps = np.array(x_values), np.array(flame_temps)
        flame_table = np.stack((x_values, flame_temps))
        return flame_table