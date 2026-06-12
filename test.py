from domain.compounds import compounds
from domain.dissociativeflametemp import DissociativeReaction
from domain.adiabatic_flame_temp import CombustionReaction
from numpy import array, float64


# test_reaction = DissociativeReaction(
#         fuels = {"Methane": 1},
#         oxi = {"Oxygen": 1},
#         temps = {"Methane": 298.15, "Oxygen": 298.15},
#         pres_bar = 1e5,
#         conc_res = 5
#     )

# # print(test_reaction._equilibrium_residual("Methane", guess, species_indices)) # This does not agree with manual calculation.)

# # methane = compounds["Methane"]
# # carbon_monoxide = compounds["Carbon_Monoxide"]
# # oxygen = compounds["Oxygen"]
# # hydrogen_monatomic = compounds["Hydrogen_Monatomic"]
# # for species in [methane, carbon_monoxide, oxygen, hydrogen_monatomic]:
# #     for temp in [300, 1000, 2000, 3000]:
# #         print(f"{species.name} logKf at {temp} K: {species.logKf(temp)}")



# guess = array([-20, -0.4771, -20, -20, -20, -20, -20, -20, -20, -0.1761, 3000.0], dtype=float64)

# # test_reaction._residual_list(guess, 1)

# print(test_reaction.equilibrate(1))

test = CombustionReaction(
    fuels = {"Methane": 1},
    oxidants = {"Oxygen": 1},
    temps = {"Methane": 298.15, "Oxygen": 298.15},
    pressure = 1e5
)

print(test._independents)
print(test._indep_matr)
print(test._item_indices)