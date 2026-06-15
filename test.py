from domain.compounds import compounds
from domain.dissociativeflametemp import DissociativeReaction
from domain.adiabatic_flame_temp import CombustionReaction
from math import log10
import numpy as np


test_reaction = CombustionReaction(
    fuels = {"Methane": 1},
    oxidants = {"Oxygen": 1},
    temps = {"Methane": 298.15, "Oxygen": 298.15},
    pressure = 1e5
)

init_atoms = np.array([1.32, 0.33, 1.34])
init_guess = np.full(8, -0.3)
init_guess[7] = 3000.0

print(test_reaction._calc_guess_vector(init_atoms, init_guess))
print(test_reaction._item_indices)
print(test_reaction._residual_indices)
print(test_reaction._residual_function(32, init_guess))

# print(np.any(init_guess <= 0))