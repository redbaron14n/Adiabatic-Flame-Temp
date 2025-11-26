# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Compound Class File
# ###################

import numpy as np
from numpy.typing import NDArray
import pandas as pd
from scipy.interpolate import make_interp_spline

DATA_FILE = "thermochemical_data.csv"
TD = pd.read_csv(DATA_FILE)

class Compound:

    ##### Initialization Methods #####

    def __init__(self, name: str, formula: str, id: str, ref_temp: float = 298.15):

        self._set_name(name)
        self._set_formula(formula)
        self._set_id(id)
        self._set_data()
        self._set_ref_temp(ref_temp)

    def _set_name(self, name: str):

        self.name = name

    def _set_formula(self, formula: str):

        self.formula = formula

    def _set_id(self, id: str):

        self._id = id

    def _set_data(self):

        self.data_table = TD[TD["Compound"] == self._id]
        self._set_temperatures(self.data_table["T"])
        self._set_SH(self.data_table["SH"])
        self._set_Hf(self.data_table["Hf"])
        self._set_logKf(self.data_table["logKf"])

    def _set_temperatures(self, list: pd.Series):

        temp_array = list.to_numpy(dtype=np.float64)
        self._temperatures = temp_array

    def _set_SH(self, list: pd.Series):

        SH_array = list.to_numpy()
        self._SH_list = SH_array
        self._SH_function = make_interp_spline(self.get_temperatures(), SH_array, k=1)

    def _set_Hf(self, list: pd.Series):

        Hf_array = list.to_numpy()
        self._Hf_list = Hf_array
        self._Hf_function = make_interp_spline(self.get_temperatures(), Hf_array, k=1)

    def _set_logKf(self, list: pd.Series):
        
        logKf_array = list.to_numpy()
        logKf_array, logKf_finite_array = self._get_numerical(logKf_array)
        self._logKf_list = logKf_array
        self._logKf_function = make_interp_spline(self.get_temperatures(), logKf_finite_array, k=1)

    def _set_ref_temp(self, ref_temp: float):

        self.ref_temp = ref_temp
        self.stdHf = self.Hf(ref_temp)

    ##### Attribute Methods #####

    def get_temperatures(self) -> NDArray[np.float64]:

        temp_array = self._temperatures
        return temp_array
    
    def get_SH_list(self) -> NDArray[np.float64]:

        SH_array = self._SH_list
        return SH_array
    
    def get_Hf_list(self) -> NDArray[np.float64]:

        Hf_array = self._Hf_list
        return Hf_array
    
    def get_logKf_list(self) -> NDArray[np.float64]:

        logKf_list = self._logKf_list
        return logKf_list
    
    ##### Functional Methods #####

    def SH(self, temperature: float) -> float:

        """
        Returns the sensible heat (kJ/mol) of the compound at a given temperature (K).
        """

        value = float(self._SH_function(temperature))
        return value
    
    def Hf(self, temperature: float) -> float:

        """
        Returns the heat of formation (kJ/mol) of the compound at a given temperature (K).
        """

        value = float(self._Hf_function(temperature))
        return value
    
    def logKf(self, temperature: float) -> float:

        """
        Returns the logKf of the compound at a given temperature (K).
        """

        value = self._logKf_function(temperature)
        return value
    
    def _get_finite_list(self, list: NDArray) -> NDArray:

        finite_list = np.copy(list)
        finite_mask = np.isfinite(list)
        max_finite = np.max(list[finite_mask])
        finite_list[~finite_mask] = max_finite * 1e6
        return finite_list

    def _get_numerical(self, list: NDArray) -> tuple[NDArray, NDArray]:

        numerical_list = np.copy(list)
        for i in range(len(list)):
            if list[i] == "inf":
                numerical_list[i] = np.inf
        finite_list = self._get_finite_list(numerical_list)
        return numerical_list, finite_list