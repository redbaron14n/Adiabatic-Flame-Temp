# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# Compound Class File
# ###################

import numpy as np
from numpy.typing import NDArray
from domain.compound_data import CompoundData
from scipy.interpolate import BSpline, make_interp_spline

STANDARD_REF_TEMP = 298.15


class Compound:

    def __init__(self, name: str, formula: str, id: str, data: CompoundData):
        self.name: str = name
        self.formula: str = formula
        self.id: str = id
        self._data: CompoundData = data

        self._Cf_function = make_interp_spline(
            self._data.temperatures,
            self._data.Cf_list,
            k=1,
        )

        self._S_function = make_interp_spline(
            self._data.temperatures,
            self._data.S_list,
            k=1,
        )

        self._DS_function = self._make_finite_function(self._data.DS_list)

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

        self._Gf_function = make_interp_spline(
            self._data.temperatures,
            self._data.Gf_list,
            k=1,
        )

        self._logKf_function = self._make_finite_function(self._data.logKf_list)

        self.stdHf = self._Hf_function(STANDARD_REF_TEMP)


    def Cf(self, temperature: float) -> float:

        """
        Returns the heat capacity (kJ/mol-K) of the compound at a given temperature (K).
        """

        value = float(self._Cf_function(temperature))
        return value
    

    def S(self, temperature: float) -> float:

        """
        Returns the entropy (kJ/mol-K) of the compound at a given temperature (K).
        """

        value = float(self._S_function(temperature))
        return value
    

    def DS(self, temperature: float) -> float:

        """
        Returns the change in entropy (kJ/mol-K) of the compound at a given temperature (K).
        """

        value = float(self._DS_function(temperature))
        return value


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


    def Gf(self, temperature: float) -> float:

        """
        Returns the Gibbs free energy of formation (kJ/mol) of the compound at a given temperature (K).
        """

        value = float(self._Gf_function(temperature))
        return value


    def logKf(self, temperature: float) -> float:

        """
        Returns the logKf of the compound at a given temperature (K).
        """

        value = self._logKf_function(temperature)
        return value


    def get_temperatures(self) -> NDArray:

        """
        Returns the list of temperatures (K) for which data is available.
        """

        return self._data.temperatures


    def get_data(self, label: str):

        """
        Returns the data list corresponding to the given label.
        Labels: "Cf", "S", "DS", "Hf", "SH", "Gf", "logKf"
        """

        match label:
            case "Cf":
                return self._data.Cf_list
            case "S":
                return self._data.S_list
            case "DS":
                return self._data.DS_list
            case "Hf":
                return self._data.Hf_list
            case "SH":
                return self._data.SH_list
            case "Gf":
                return self._data.Gf_list
            case "logKf":
                return self._data.logKf_list
            case _:
                raise ValueError(f"Data label '{label}' not recognized.")

    def _make_finite_function(self, list: NDArray) -> BSpline:

        """
        Separate function maker method as DS and logKf tables include np.inf values.
        """

        finite_list: NDArray = self._get_finite_list(list)
        return make_interp_spline(
            self._data.temperatures,
            finite_list,
            k=1,
        )

    def _get_finite_list(self, list: NDArray) -> NDArray:

        """
        Turns np.inf values into 1e6 * largest finite value
        """

        finite_list = np.copy(list)
        finite_mask = np.isfinite(list)
        max_finite = np.max(list[finite_mask])
        finite_list[~finite_mask] = max_finite * 1e6
        return finite_list
