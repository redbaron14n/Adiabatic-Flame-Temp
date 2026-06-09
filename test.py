from domain.compounds import compounds
from domain.dissociativeflametemp import DissociativeReaction
from numpy import array


guess = array([3000.0, 0.5, 0.5, 0.1, 0.4, 0.7, 0.2, 0.9, 0.1, 0.5, 0.1])
species_indices = {"T": 0, "Methane": 1, "Oxygen": 2, "Carbon_Dioxide": 3, "Water": 4, "Carbon_Monoxide": 5, "Oxygen_Monatomic": 6, "Carbon": 7, "Hydrogen": 8, "Hydroxyl": 9, "Hydrogen_Monatomic": 10}

guess_actual = array([0.9, 0.1, 0.7, 0.1, 0.1, 0.5, 0.5, 0.5, 0.2, 0.4, 3000.0])


test_reaction = DissociativeReaction(
        fuels = {"Methane": 1},
        oxi = {"Oxygen": 1},
        temps = {"Methane": 300, "Oxygen": 300},
        pres_bar = 1.0,
        conc_res = 5
    )

# print(test_reaction._equilibrium_residual("Methane", guess, species_indices)) # This does not agree with manual calculation.)

print(test_reaction.equilibrate(1))

# methane = compounds["Methane"]
# carbon_monoxide = compounds["Carbon_Monoxide"]
# oxygen = compounds["Oxygen"]
# hydrogen_monatomic = compounds["Hydrogen_Monatomic"]
# for species in [methane, carbon_monoxide, oxygen, hydrogen_monatomic]:
#     for temp in [300, 1000, 2000, 3000]:
#         print(f"{species.name} logKf at {temp} K: {species.logKf(temp)}")

for diss_obj in test_reaction._dissociations.values():
    print(f"{diss_obj.molecule_id} dissociates into {diss_obj._stoich}")