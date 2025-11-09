# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Thermochemical Data File Editor
# ###################

import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import make_interp_spline

TEMPERATURE_LIST = np.array([0, 100, 200, 250, 298.15, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800, 5900, 6000])

class Compound:

    def __init__(self, formula: str, sensible_heat: NDArray[np.float64], heat_of_formation: NDArray[np.float64], logKf: NDArray[np.float64], temperature: NDArray[np.float64] = TEMPERATURE_LIST):

        self.__set_formula(formula)
        self.__set_temp_list(temperature)
        self.__set_sh_list(sensible_heat)
        self.__set_hf_list(heat_of_formation)
        self.__set_logKf_list(logKf)

    def __set_formula(self, formula: str):

        self.__formula = formula

    def formula(self) -> str:

        formula = self.__formula
        return formula

    def __set_temp_list(self,t: NDArray[np.float64]):

        if len(np.shape(t)) != 1:
            raise ValueError("Temperature array should be one-dimensional.")
        self.__temp_list = t

    def __get_temp_list(self) -> NDArray[np.float64]:

        temp = self.__temp_list
        return temp

    def __set_sh_list(self, specific_heat: NDArray[np.float64]):

        if np.shape(specific_heat) != np.shape(self.__get_temp_list()):
            raise ValueError("List of specific heats has differing length than temperature list.")
        self.__sh_list = specific_heat
        self.__sh_function = make_interp_spline(self.__get_temp_list(), specific_heat, k=1)

    def __get_sh_list(self) -> NDArray[np.float64]:

        sh = self.__sh_list
        return sh
    
    def sh(self, temperature: float) -> float:

        sh = self.__sh_function(temperature)
        return sh

    def __set_hf_list(self, hf: NDArray[np.float64]):

        if np.shape(hf) != np.shape(self.__get_temp_list()):
            raise ValueError("List of heats of formation has differing length than temperature list.")
        self.__hf_list = hf
        self.__hf_function = make_interp_spline(self.__get_temp_list(), hf, k=1)

    def __get_hf_list(self) -> NDArray[np.float64]:

        hf = self.__hf_list
        return hf
    
    def hf(self, temperature: float) -> float:

        hf = self.__hf_function(temperature)
        return hf
        
    def __set_logKf_list(self, logKf: NDArray[np.float64]):

        if np.shape(logKf) != np.shape(self.__get_temp_list()):
            raise ValueError("List of log Kf values has differing length than temperature list.")
        self.__logKf_list = logKf
        # self.__logKf_function = make_interp_spline(self.__get_temp_list(), logKf, k=1) # ADDRESS INF VALUE

    def __get_logKf_list(self) -> NDArray[np.float64]:

        logKf = self.__logKf_list
        return logKf

    # def logKf(self, temperature: float) -> float:

    #     logKf = self.__logKf_function(temperature)
    #     return logKf

    def hf_table(self, transpose: bool = False) -> NDArray[np.float64]:

        table = np.stack((self.__get_temp_list(), self.__get_hf_list()))
        if transpose:
            table = np.transpose(table)
        return table
    
    def sh_table(self, transpose: bool = False) -> NDArray[np.float64]: # Output is formatted odd. Investigate

        table = np.stack((self.__get_temp_list(), self.__get_sh_list()))
        if transpose:
            table = np.transpose(table)
        return table
    
    def logKf_table(self, transpose: bool = False) -> NDArray[np.float64]:

        table = np.stack((self.__get_temp_list(), self.__get_logKf_list()))
        if transpose:
            table = np.transpose(table)
        return table

##### Methane (CH4) #####

methane_sh = np.array([-10.024, -6.698, -3.368, -1.679, 0, 0.066, 1.903, 3.861, 5.957, 8.200, 13.130, 18.635, 24.675, 31.205, 38.179, 45.549, 53.270, 61.302, 69.608, 78.153, 86.910, 95.853, 104.960, 114.212, 123.592, 133.087, 142.684, 152.371, 162.141, 171.684, 181.893, 191.862, 201.885, 211.958, 222.076, 232.235, 242.431, 252.662, 262.925, 273.217, 283.536, 293.881, 304.248, 314.637, 325.045, 335.473, 345.918, 356.379, 366.855, 377.345, 387.849, 398.365, 408.893, 419.432, 429.982, 440.541, 451.110, 461.688, 472.274, 482.867, 493.469, 504.077, 514.692, 525.314, 535.942])
methane_hf = np.array([-66.911, -69.644, -72.027, -73.426, -74.873, -74.929, -76.461, -77.969, -79.422, -80.802, -83.308, -85.452, -87.238, -88.692, -89.849, -90.750, -91.437, -91.945, -92.308, -92.553, -92.703, -92.780, -92.797, -92.770, -92.709, -92.624, -92.521, -92.409, -92.291, -92.174, -92.060, -91.954, -91.857, -91.773, -91.705, -91.653, -91.621, -91.609, -91.619, -91.654, -91.713, -91.798, -91.911, -92.051, -92.222, -92.422, -92.652, -92.914, -93.208, -93.533, -93.891, -94.281, -94.702, -95.156, -95.641, -96.157, -96.703, -97.278, -97.882, -98.513, -99.170, -99.852, -100.557, -101.284, -102.032])
methane_logkf = np.array([np.inf, 33.615, 15.190, 11.395, 8.894, 8.813, 6.932, 5.492, 4.350, 3.420, 1.993, 0.943, 0.138, -0.500, -1.018, -1.447, -1.807, -2.115, -2.379, -2.609, -2.810, -2.989, -3.147, -3.289, -3.416, -3.531, -3.636, -3.732, -3.819, -3.899, -3.973, -4.042, -4.105, -4.164, -4.219, -4.271, -4.319, -4.364, -4.407, -4.447, -4.485, -4.521, -4.555, -4.588, -4.619, -4.648, -4.676, -4.703, -4.729, -4.753, -4.777, -4.800, -4.822, -4.843, -4.863, -4.883, -4.902, -4.920, -4.938, -4.955, -4.972, -4.988, -5.004, -5.019, -5.034])

methane = Compound("CH4", methane_sh, methane_hf, methane_logkf)