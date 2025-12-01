from chempy import balance_stoichiometry

reactants = {"CH4", "O2", "N2", "Ar"}
products = {"CO2", "H2O"}
balanced_reactants, balanced_products = balance_stoichiometry(reactants, products)
print(balanced_reactants, balanced_products)