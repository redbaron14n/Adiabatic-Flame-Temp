import numpy as np
import pandas as pd

from numpy.typing import NDArray
from scipy.interpolate import make_interp_spline

from domain.compound_data import CompoundData

STANDARD_REF_TEMP = 298.15

DATA_FILE = "thermochemical_data.csv"
TD: pd.DataFrame = pd.read_csv(DATA_FILE)


class CompoundLoader:
    def load(self, id) -> CompoundData:
        data_table = TD[TD["Compound"] == id]

        temperatures = data_table["T"].to_numpy(dtype=np.float64)
        SH_list = data_table["SH"].to_numpy()
        Hf_list = data_table["Hf"].to_numpy()
        logKf_list = self._convert_infs(data_table["logKf"])

        compound_data = CompoundData(
            temperatures=temperatures,
            SH_list=SH_list,
            Hf_list=Hf_list,
            logKf_list=logKf_list,
        )

        return compound_data

    def _convert_infs(self, series: NDArray) -> NDArray:
        converted_list = np.copy(series)
        for i in range(len(series)):
            if series[i] == "inf":
                converted_list[i] = np.inf
        return converted_list
