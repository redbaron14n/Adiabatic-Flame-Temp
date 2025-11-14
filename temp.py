# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Thermochemical Data File Editor
# ###################

import numpy as np
from numpy.typing import NDArray
import pandas as pd
from scipy.interpolate import make_interp_spline

DATA_FILE = "thermochemical_data.csv"
TD = pd.read_csv(DATA_FILE)

def get_finite_list(list: NDArray) -> NDArray:

    finite_list = np.copy(list)
    finite_mask = np.isfinite(list)
    max_finite = np.max(list[finite_mask])
    finite_list[~finite_mask] = max_finite * 1e6
    return finite_list

def get_numerical(list: NDArray) -> tuple[NDArray, NDArray]:

    numerical_list = np.copy(list)
    for i in range(len(list)):
        if list[i] == "inf":
            numerical_list[i] = np.inf
    finite_list = get_finite_list(numerical_list)
    return numerical_list, finite_list

class Compound:

    ##### Initialization Methods #####

    def __init__(self, name: str, formula: str, id: str):

        self.__set_name(name)
        self.__set_formula(formula)
        self.__set_id(id)
        self.__set_data(id)

    def __set_name(self, name: str):

        self.name = name

    def __set_formula(self, formula: str):

        self.formula = formula

    def __set_id(self, id: str):

        self.__id = id

    def __set_data(self, id: str):

        self.data_table = TD[TD["Compound"] == self.__id]
        self.__set_temperatures(self.data_table["T"])
        self.__set_SH(self.data_table["SH"])
        self.__set_Hf(self.data_table["Hf"])
        self.__set_logKf(self.data_table["logKf"])

    def __set_temperatures(self, list: pd.Series):

        temp_array = list.to_numpy(dtype=np.float64)
        self.__temperatures = temp_array

    def __set_SH(self, list: pd.Series):

        SH_array = list.to_numpy()
        self.__SH_list = SH_array
        self.__SH_function = make_interp_spline(self.get_temperatures(), SH_array, k=1)

    def __set_Hf(self, list: pd.Series):

        Hf_array = list.to_numpy()
        self.__Hf_list = Hf_array
        self.__Hf_function = make_interp_spline(self.get_temperatures(), Hf_array, k=1)

    def __set_logKf(self, list: pd.Series):
        
        logKf_array = list.to_numpy()
        logKf_array, logKf_finite_array = get_numerical(logKf_array)
        self.__logKf_list = logKf_array
        self.__logKf_function = make_interp_spline(self.get_temperatures(), logKf_finite_array, k=1)

    ##### Attribute Methods #####

    def get_temperatures(self) -> NDArray[np.float64]:

        temp_array = self.__temperatures
        return temp_array
    
    def get_SH_list(self) -> NDArray[np.float64]:

        SH_array = self.__SH_list
        return SH_array
    
    def get_Hf_list(self) -> NDArray[np.float64]:

        Hf_array = self.__Hf_list
        return Hf_array
    
    def get_logKf_list(self) -> NDArray[np.float64]:

        logKf_list = self.__logKf_list
        return logKf_list
    
    ##### Functional Methods #####

    def SH(self, temperature: float) -> float:

        """
        Returns the sensible heat (kJ/mol) of the compound at a given temperature (K).
        """

        value = self.__SH_function(temperature)
        return value
    
    def Hf(self, temperature: float) -> float:

        """
        Returns the heat of formation (kJ/mol) of the compound at a given temperature (K).
        """

        value = self.__Hf_function(temperature)
        return value
    
    def logKf(self, temperature: float) -> float:

        """
        Returns the logKf of the compound at a given temperature (K).
        """

        value = self.__logKf_function(temperature)
        return value

Methane = Compound("Methane", "CH4", "Methane")