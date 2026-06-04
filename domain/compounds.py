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

compounds["Oxygen_Monatomic"] = Compound(
    name = "Monatomic Oxygen",
    formula = "O",
    id = "Oxygen_Monatomic",
    data = load_compound_data("Oxygen_Monatomic")
)

compounds["Carbon"] = Compound(
    name = "Carbon",
    formula = "C",
    id = "Carbon",
    data = load_compound_data("Carbon"),
    state = "s"
)

compounds["Hydrogen_Monatomic"] = Compound(
    name = "Monatomic Hydrogen",
    formula = "H",
    id = "Hydrogen_Monatomic",
    data = load_compound_data("Hydrogen_Monatomic")
)

compounds["Nitrogen_Oxide"] = Compound(
    name = "Nitrogen Oxide",
    formula = "NO",
    id = "Nitrogen_Oxide",
    data = load_compound_data("Nitrogen_Oxide")
)

compounds["Hydroxyl"] = Compound(
    name = "Hydroxyl",
    formula = "OH",
    id = "Hydroxyl",
    data = load_compound_data("Hydroxyl")
)

compounds["Carbon_Monoxide"] = Compound(
    name = "Carbon Monoxide",
    formula = "CO",
    id = "Carbon_Monoxide",
    data = load_compound_data("Carbon_Monoxide")
)