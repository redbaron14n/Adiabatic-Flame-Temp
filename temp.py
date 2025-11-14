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

DATA_FILE = "thermochemical_data.csv"
TD = pd.read_csv(DATA_FILE)

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
        self.__set_logKf(self.data_table["logKf"])

    def __set_temperatures(self, list: pd.Series):

        temp_array = list.to_numpy(dtype=np.float64)
        self.__temperatures = temp_array

    def __set_logKf(self, list: pd.Series):
        
        logKf_array = list.to_numpy()
        logKF_array, logKf_finite_array = get_numerical(logKF_array)
        self.__logKf_list = logKF_array
        self.logKf_function() # WIP

    ##### Attribute Methods #####



    def get_temperatures(self) -> NDArray[np.float64]:

        temp_array = self.__temperatures
        return temp_array

Methane = Compound("Methane", "CH4", "Methane")
print(Methane.data_table)