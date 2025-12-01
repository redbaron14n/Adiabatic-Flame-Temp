# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Configuration File
# ###################


def products_from_reactants(reactants: set[str], dissociation: bool = False) -> set[str]:  # Placeholder function for potential reaction product generation

    if (reactants == {"Methane", "Oxygen"}) and not dissociation:
        products = {"Carbon_Dioxide", "Water"}
    elif (reactants == {"Hydrogen", "Oxygen"}) and not dissociation:
        products = {"Water"}
    else:
        raise NotImplementedError(
            f"Reaction product generation not implemented for given reactants.\n"
            f"{reactants}"
        )
    return products