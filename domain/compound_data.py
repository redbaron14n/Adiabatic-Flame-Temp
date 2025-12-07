import numpy
from dataclasses import dataclass

# Dataclass stores data pulled from .csv file as named attributes

@dataclass
class CompoundData:
    temperatures: numpy.ndarray
    SH_list: numpy.ndarray
    Hf_list: numpy.ndarray
    logKf_list: numpy.ndarray
