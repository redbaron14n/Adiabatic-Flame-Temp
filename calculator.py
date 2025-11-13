# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Calculator File
# ###################

from ThermochemicalDataEditor import *

##### Reaction Determination #####

# Replace with an actual method of doing this lmfao

def determine_products(reactants: tuple[Compound]) -> tuple[Compound]:

    if reactants == (Methane, Oxygen):
        products = (Carbon_Dioxide, Water) # Need to actually define this in the data still
    return products