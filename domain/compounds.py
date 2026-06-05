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
    dissociates = {"Carbon_Monoxide", "Oxygen_Monatomic", "Carbon", "Oxygen"},
    composition = {"Carbon": 1, "Oxygen": 1}
)

compounds["Methane"] = Compound(
    name="Methane",
    formula="CH4",
    id="Methane",
    data=load_compound_data("Methane"),
    composition = {"Carbon": 1, "Hydrogen": 2}
)

compounds["Water"] = Compound(
    name="Water",
    formula="H2O",
    id="Water",
    data=load_compound_data("Water"),
    dissociates = {"Hydrogen", "Oxygen_Monatomic", "Hydroxyl", "Hydrogen_Monatomic"},
    composition = {"Hydrogen": 1, "Oxygen": 0.5}
)

compounds["Oxygen"] = Compound(
    name="Oxygen",
    formula="O2",
    id="Oxygen",
    data=load_compound_data("Oxygen"),
    dissociates = {"Oxygen_Monatomic"}
)

compounds["Hydrogen"] = Compound(
    name = "Hydrogen",
    formula = "H2",
    id = "Hydrogen",
    data = load_compound_data("Hydrogen"),
    dissociates = {"Hydrogen_Monatomic"}
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
    data = load_compound_data("Oxygen_Monatomic"),
    composition = {"Oxygen": 0.5}
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
    data = load_compound_data("Hydrogen_Monatomic"),
    composition = {"Hydrogen": 0.5}
)

compounds["Nitrogen_Oxide"] = Compound(
    name = "Nitrogen Oxide",
    formula = "NO",
    id = "Nitrogen_Oxide",
    data = load_compound_data("Nitrogen_Oxide"),
    composition = {"Nitrogen": 0.5, "Oxygen": 0.5}
)

compounds["Hydroxyl"] = Compound(
    name = "Hydroxyl",
    formula = "OH",
    id = "Hydroxyl",
    data = load_compound_data("Hydroxyl"),
    dissociates = {"Oxygen_Monatomic", "Hydrogen_Monatomic"},
    composition = {"Hydrogen": 0.5, "Oxygen": 0.5}
)

compounds["Carbon_Monoxide"] = Compound(
    name = "Carbon Monoxide",
    formula = "CO",
    id = "Carbon_Monoxide",
    data = load_compound_data("Carbon_Monoxide"),
    dissociates = {"Carbon", "Oxygen_Monatomic"},
    composition = {"Carbon": 1, "Oxygen": 0.5}
)

compounds_by_formula: dict[str, Compound] = {c.formula: c for c in compounds.values()}