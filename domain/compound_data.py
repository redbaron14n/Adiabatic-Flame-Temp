# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# CompoundData Dataclass File
# ###################

import numpy
from dataclasses import dataclass


"""
Dataclass stores data pulled from .csv file as named attributes
"""
@dataclass
class CompoundData:
    temperatures: numpy.ndarray
    Cf_list: numpy.ndarray
    S_list: numpy.ndarray
    DS_list: numpy.ndarray
    SH_list: numpy.ndarray
    Hf_list: numpy.ndarray
    Gf_list: numpy.ndarray
    logKf_list: numpy.ndarray
