import numpy
import services.comp_loader as loader
from dataclasses import dataclass


@dataclass
class CompoundData:
    temperatures: numpy.ndarray
    SH_list: numpy.ndarray
    Hf_list: numpy.ndarray
    logKf_list: numpy.ndarray
