from services.comp_loader import CompoundLoader
from domain.compound_data import CompoundData
from domain.compound2 import Compound2


def load_compound_data(compound_id: str) -> CompoundData:
    # replace CompoundLoader with other types of loaders as needed
    loader: CompoundLoader = CompoundLoader()
    return loader.load(compound_id)  # eg "Carbon_Dioxide"


compounds: dict = {}

compounds["Carbon_Dioxide"] = Compound2(
    name="Carbon Dioxide",
    formula="CO2",
    id="Carbon_Dioxide",
    data=load_compound_data("Carbon_Dioxide"),
)

compounds["Methane"] = Compound2(
    name="Methane",
    formula="CH4",
    id="Methane",
    data=load_compound_data("Methane"),
)

compounds["Water"] = Compound2(
    name="Water",
    formula="H2O",
    id="Water",
    data=load_compound_data("Water"),
)

compounds["Oxygen"] = Compound2(
    name="Oxygen",
    formula="O2",
    id="Oxygen",
    data=load_compound_data("Oxygen"),
)
