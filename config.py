# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Configuration File
# ###################


def products_from_reactants(reactants: set[str], dissociation: bool = False) -> tuple[set[str], set[str]]:  # Placeholder function for potential reaction product generation

    inert_reactants = set()
    if (reactants == {"Methane", "Oxygen"}) and not dissociation:
        products = {"Carbon_Dioxide", "Water"}
    elif (reactants == {"Hydrogen", "Oxygen"}) and not dissociation:
        products = {"Water"}
    elif (reactants == {"Methane", "Nitrogen", "Oxygen", "Argon"}) and not dissociation:
        products = {"Carbon_Dioxide", "Water"}
        inert_reactants = {"Nitrogen", "Argon"}
    elif (reactants == {"Methane", "Nitrogen", "Oxygen"}) and not dissociation:
        products = {"Carbon_Dioxide", "Water"}
        inert_reactants = {"Nitrogen"}
    elif (reactants == {"Methane", "Oxygen", "Argon"}) and not dissociation:
        products = {"Carbon_Dioxide", "Water"}
        inert_reactants = {"Argon"}
    elif (reactants == {"Hydrogen", "Nitrogen", "Oxygen", "Argon"}) and not dissociation:
        products = {"Water"}
        inert_reactants = {"Nitrogen", "Argon"}
    elif (reactants == {"Hydrogen", "Nitrogen", "Oxygen"}) and not dissociation:
        products = {"Water"}
        inert_reactants = {"Nitrogen"}
    elif (reactants == {"Hydrogen", "Oxygen", "Argon"}) and not dissociation:
        products = {"Water"}
        inert_reactants = {"Argon"}
    else:
        raise NotImplementedError(
            f"Reaction product generation not implemented for given reactants.\n"
            f"{reactants}"
        )
    return products, inert_reactants