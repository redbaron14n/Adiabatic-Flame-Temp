# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Configuration File
# ###################

from class_files.compound import Compound

Compound_list = []
CarbonDioxide = Compound("Carbon Dioxide", "CO2", "Carbon_Dioxide")
Compound_list.append(CarbonDioxide)
Methane = Compound("Methane", "CH4", "Methane")
Compound_list.append(Methane)
Oxygen = Compound("Oxygen", "O2", "Oxygen")
Compound_list.append(Oxygen)
Water = Compound("Water", "H2O", "Water")
Compound_list.append(Water)

def products_from_reactants(reactants: set[Compound], dissociation: bool) -> set[Compound]: # Placeholder function for potential reaction product generation

    if (reactants == {Methane, Oxygen}) and not dissociation:
        products = {CarbonDioxide, Water}
    else:
        raise NotImplementedError("Reaction product generation not implemented for given reactants.")
    return products