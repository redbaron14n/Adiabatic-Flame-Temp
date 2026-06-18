from domain.compounds import compounds
from domain.dissociation import Dissociation
from domain.dissociativeflametemp import DissociativeReaction
from domain.adiabatic_flame_temp import CombustionReaction
from domain.reaction import Reaction
from math import log10
import numpy as np



test_reaction = CombustionReaction(
    fuels = {"Hydrogen": 1},
    oxidants = {"Oxygen": 1},
    temps = {"Hydrogen": 298.15, "Oxygen": 298.15},
    pressure = 1e5,
    conc_res = 5
)

test_basic_reaction = Reaction(
    reactants = {"Hydrogen", "Oxygen"},
    temperatures = {"Hydrogen": 298.15, "Oxygen": 298.15}
)

# print(test_reaction._conc_list)
# print(test_reaction._dependents)
# print(test_reaction._item_indices)
# print(test_reaction._present_atoms)
# print(test_reaction._independents)
# print(test_reaction._indep_matr)
# print(test_reaction._residual_indices)
# print(test_reaction._dissociations)
# print(test_reaction._stoich)
print(test_reaction._reactants)
# print(test_reaction._basic_products)

test_conc = test_reaction._conc_list[3]
print(test_basic_reaction.calc_flame_table(test_reaction._conc_list))