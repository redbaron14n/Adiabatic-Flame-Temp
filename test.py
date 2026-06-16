from domain.compounds import compounds
from domain.dissociation import Dissociation
from domain.dissociativeflametemp import DissociativeReaction
from domain.adiabatic_flame_temp import CombustionReaction
from math import log10
import numpy as np


test_reaction = CombustionReaction(
    fuels = {"Methane": 1},
    oxidants = {"Oxygen": 1},
    temps = {"Methane": 298.15, "Oxygen": 298.15},
    pressure = 1e5,
    conc_res = 99
)

init_atoms = np.array([1.32, 0.33, 1.34])

init_log_guess = test_reaction._calc_init_log_guess(test_reaction._conc_list[32])
print(init_log_guess, init_log_guess.shape)
bounds = test_reaction._bounds
print(bounds, bounds[0].shape, bounds[1].shape)
full_real_guess = test_reaction._calc_guess_vector(init_atoms, init_log_guess)
print(full_real_guess, full_real_guess.shape)
log_guess = test_reaction._convert_guess_to_log(full_real_guess)
print(log_guess, log_guess.shape)
print(test_reaction._item_indices)
print(test_reaction._dissociations)