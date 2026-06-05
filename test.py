from domain.dissociativeflametemp import DissociativeReaction
from numpy import array


guess = array([750.0, 0.5, 0.5, 0.1, 0.4, 0.7, 0.2, 0.9, 0.1, 0.5, 0.1])
species_indices = {"T": 0, "Methane": 1, "Oxygen": 2, "Carbon_Dioxide": 3, "Water": 4, "Carbon_Monoxide": 5, "Oxygen_Monatomic": 6, "Carbon": 7, "Hydrogen": 8, "Hydroxyl": 9, "Hydrogen_Monatomic": 10}


test_reaction = DissociativeReaction(
        fuels = {"Methane": 1},
        oxi = {"Oxygen": 1},
        temps = {"Methane": 300, "Oxygen": 300},
        pres_bar = 1.0,
        conc_res = 5
    )

print(test_reaction._equilibrium_residual("Methane", guess, species_indices)) # This does not agree with manual calculation.