# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Configuration File
# ###################

INERTS = {"Argon", "Nitrogen"}

def products_from_reactants(reactants: set[str], dissociation: bool = False) -> tuple[set[str], set[str]]:  # Placeholder function for potential reaction product generation

    active_reactants = reactants - INERTS
    inert_reactants = reactants & INERTS
    if (active_reactants == {"Methane", "Oxygen"}) and not dissociation:
        products = {"Carbon_Dioxide", "Water"}
    elif (active_reactants == {"Hydrogen", "Oxygen"}) and not dissociation:
        products = {"Water"}
    else:
        raise NotImplementedError(
            f"Reaction product generation not implemented for given reactants.\n"
            f"{reactants}"
        )
    return products, inert_reactants