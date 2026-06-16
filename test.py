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

conc_dict = test_reaction._conc_list[32]
init_atoms = np.array([1.32, 0.33, 1.34])

# init_log_guess = test_reaction._calc_init_log_guess(conc_dict)
# print(f"Initial log guess: {init_log_guess}, shape: {init_log_guess.shape}")
# bounds = test_reaction._bounds
# print(f"Bounds: {bounds}, lower shape: {bounds[0].shape}, upper shape: {bounds[1].shape}")
# full_real_guess = test_reaction._calc_guess_vector(init_atoms, init_log_guess)
# print(f"Full real guess: {full_real_guess}, shape: {full_real_guess.shape}")
# log_guess = test_reaction._convert_guess_to_log(full_real_guess)
# print(f"Log guess: {log_guess}, shape: {log_guess.shape}")
# print(f"Item indices: {test_reaction._item_indices}")
# print(f"Residual indices: {test_reaction._residual_indices}")
# print(test_reaction.equilibrate(conc_dict))

# real_init_guess = np.array([1.00000e-50, 3.31517e-01, 3.11992e-03, 9.04888e-20, 3.02218e-13, 1.31646e-12, 5.67045e-05, 5.61426e+03, -4.63654e-03, 2.25853e-03, 6.69273e-01], dtype=np.float64)
# log_input1 = np.array([-50, -4.79494e-01, -2.50586, -19.04341, -12.51968, -11.88059, -4.24638, 5.61426e+03], dtype=np.float64)
# log_input2 = np.copy(log_input1)
# delta_co = 1e-3
# log_input2[2] += delta_co
# residual1 = test_reaction._residual_function(log_input1, init_atoms, conc_dict)
# residual2 = test_reaction._residual_function(log_input2, init_atoms, conc_dict)
# print(residual1)
# print(residual2)
# print((residual2 - residual1)/delta_co)

log_guess = np.array([-50, -5.37101e-01, -1.40159, -2.03070e+01, -1.19567e+01, -1.10989e+01, -3.36833, 5.60963e+03], dtype=np.float64)
methane = test_reaction._calc_guess_vector(init_atoms, log_guess)[8]
print(f"Original methane: {methane}")
for i in range(len(log_guess)):
    xp = log_guess.copy()
    delta = 1e-7*(1 + abs(log_guess[i]))
    xp[i] += delta
    methane = test_reaction._calc_guess_vector(init_atoms, xp)[8]
    print(f"Perturbed {i} to {xp[i]}: {methane}")