import numpy as np
from numpy.typing import NDArray
from domain.compound_data import CompoundData
from scipy.interpolate import make_interp_spline

STANDARD_REF_TEMP = 298.15


class Compound:

    def __init__(self, name: str, formula: str, id: str, data: CompoundData):
        self.name: str = name
        self.formula: str = formula
        self.id: str = id
        self._data: CompoundData = data

        self._SH_function = make_interp_spline(
            self._data.temperatures,
            self._data.SH_list,
            k=1,
        )

        self._Hf_function = make_interp_spline(
            self._data.temperatures,
            self._data.Hf_list,
            k=1,
        )

        self._logKf_function = self._make_logKf_function()

        self.stdHf = self._Hf_function(STANDARD_REF_TEMP)

    def Hf(self, temperature: float = STANDARD_REF_TEMP) -> float:
        """
        Returns the heat of formation (kJ/mol) of the compound at a given temperature (K).
        """

        value = float(self._Hf_function(temperature))
        return value

    def SH(self, temperature: float) -> float:
        """
        Returns the sensible heat (kJ/mol) of the compound at a given temperature (K).
        """

        value = float(self._SH_function(temperature))
        return value

    def get_temperatures(self) -> NDArray:
        return self._data.temperatures

    def logKf(self, temperature: float) -> float:
        """
        Returns the logKf of the compound at a given temperature (K).
        """

        value = self._logKf_function(temperature)
        return value

    def get_data(self, label: str):
        match label:
            case "SH":
                return self._data.SH_list
            case "Hf":
                return self._data.Hf_list
            case "logKf":
                return self._data.logKf_list
            case _:
                raise ValueError(f"Data label '{label}' not recognized.")

    def _make_logKf_function(self):
        finite_logKf_list: NDArray = self._get_finite_list(self._data.logKf_list)
        return make_interp_spline(
            self._data.temperatures,
            finite_logKf_list,
            k=1,
        )

    def _get_finite_list(self, list: NDArray) -> NDArray:
        finite_list = np.copy(list)
        finite_mask = np.isfinite(list)
        max_finite = np.max(list[finite_mask])
        finite_list[~finite_mask] = max_finite * 1e6
        return finite_list
