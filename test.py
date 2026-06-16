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

# init_atoms = np.array([1.32, 0.33, 1.34])
# init_guess = np.array([-0.48149, -1, -1, -1, -1, -1, -0.17393, 3000.0])
# print(test_reaction._calc_guess_vector(init_atoms, init_guess))

conc = test_reaction._conc_list[32]
extent = test_reaction._calc_extnt_of_react(conc)
print(extent)
print(test_reaction._calc_basic_final_amouns(conc, extent))
print(test_reaction._calc_init_log_guess(conc))