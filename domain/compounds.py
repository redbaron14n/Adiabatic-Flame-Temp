# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Compound Dictionary File
# ###################

from services.comp_loader import CompoundLoader
from domain.compound_data import CompoundData
from domain.compound import Compound


def load_compound_data(compound_id: str) -> CompoundData:
    # replace CompoundLoader with other types of loaders as needed
    loader: CompoundLoader = CompoundLoader()
    return loader.load(compound_id)  # eg "Carbon_Dioxide"


compounds: dict[str, Compound] = {}

compounds["Carbon_Dioxide"] = Compound(
    name="Carbon Dioxide",
    formula="CO2",
    id="Carbon_Dioxide",
    data=load_compound_data("Carbon_Dioxide"),
)

compounds["Methane"] = Compound(
    name="Methane",
    formula="CH4",
    id="Methane",
    data=load_compound_data("Methane"),
)

compounds["Water"] = Compound(
    name="Water",
    formula="H2O",
    id="Water",
    data=load_compound_data("Water"),
)

compounds["Oxygen"] = Compound(
    name="Oxygen",
    formula="O2",
    id="Oxygen",
    data=load_compound_data("Oxygen"),
)

compounds["Hydrogen"] = Compound(
    name = "Hydrogen",
    formula = "H2",
    id = "Hydrogen",
    data = load_compound_data("Hydrogen"),
)

compounds["Nitrogen"] = Compound(
    name = "Nitrogen",
    formula = "N2",
    id = "Nitrogen",
    data = load_compound_data("Nitrogen"),
)

compounds["Argon"] = Compound(
    name = "Argon",
    formula = "Ar",
    id = "Argon",
    data = load_compound_data("Argon")
)