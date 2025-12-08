# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Compound Data Loader File
# ###################

from domain.compound_data import CompoundData
import numpy as np
from numpy.typing import NDArray
import pandas as pd

STANDARD_REF_TEMP = 298.15

DATA_FILE = "thermochemical_data.csv"
TD: pd.DataFrame = pd.read_csv(DATA_FILE)


"""
Reads thermochemical data from CSV file and loads it into CompoundData objects.
"""
class CompoundLoader:
    def load(self, id) -> CompoundData:
        data_table = TD[TD["Compound"] == id]

        temperatures = data_table["T"].to_numpy(dtype=np.float64)
        Cf_list = data_table["Cf"].to_numpy(dtype=np.float64)
        S_list = data_table["S"].to_numpy(dtype=np.float64)
        DS_list = self._convert_infs(data_table["(G-H)/T"].to_numpy()) # data includes 'inf' strings
        SH_list = data_table["SH"].to_numpy(dtype=np.float64)
        Hf_list = data_table["Hf"].to_numpy(dtype=np.float64)
        Gf_list = data_table["G"].to_numpy(dtype=np.float64)
        logKf_list = self._convert_infs(data_table["logKf"].to_numpy()) # data includes 'inf' strings

        compound_data = CompoundData(
            temperatures=temperatures,
            Cf_list=Cf_list,
            S_list=S_list,
            DS_list=DS_list,
            SH_list=SH_list,
            Hf_list=Hf_list,
            Gf_list=Gf_list,
            logKf_list=logKf_list,
        )

        return compound_data

    def _convert_infs(self, series: NDArray) -> NDArray:
        converted_list = np.copy(series)
        for i in range(len(series)):
            if series[i] == "inf":
                converted_list[i] = np.inf
        return converted_list
