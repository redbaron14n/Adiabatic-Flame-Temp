# ###################
# Ian Janes
# Professor Don Lipkin
# Adiabatic Flame Temperature
# Compound Class File
# ###################

import numpy as np
from numpy.typing import NDArray
from domain.compound_data import CompoundData
from scipy.interpolate import BSpline, make_interp_spline


class Compound:

    def __init__(
            self,
            name: str,
            formula: str,
            id: str,
            data: CompoundData,
            dissociates: set[str] = set(),
            composition: dict[str, float] = dict(),
            state: str = "g"
        ):

        """
        :param str name: The common name of the compound.
        :param str formula: The chemical formula of the compound.
        :param str id: A unique identifier for the compound.
        :param CompoundData data: A CompoundData object containing the thermodynamic data for the compound.
        :param set[str] dissociates: A set of compound IDs that this compound can dissociate into. Default is an empty set.
        :param dict[str, float] composition: A dictionary mapping the IDs of the elements in the compound to their stoichiometric coefficients for one mole of the compound. Default is an empty dictionary.
        :param str state: The reference state of the compound, either "g" for gas or "s" for solid. Default is "g".
        """

        self.name: str = name
        self.formula: str = formula
        self.id: str = id
        self.state: str = state
        self.dissociates = dissociates
        self._data: CompoundData = data
        self.composition = composition

        self._Cp_function = make_interp_spline(
            self._data.temperatures,
            self._data.Cp_list,
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

        self._set_std_temp()

        self.stdHf = self._Hf_function(self.std_temp)


    def _set_std_temp(self):

        """
        Sets the standard reference temperature for the compound based on the data provided.
        """

        temp = -1.0
        index = 0
        SH_data = self._data.SH_list
        while temp == -1.0:
            if SH_data[index] == 0.0:
                temp = self._data.temperatures[index]
            index += 1
        if temp == -1.0:
            raise ValueError("No standard reference temperature found in data.")
        self.std_temp = temp


    def Cp(self, temperature: float) -> float:

        """
        Returns the constant-pressure specific heat capacity (kJ/mol-K) of the compound at a given temperature (K).
        """

        value = float(self._Cp_function(temperature))
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
        Labels: "Cp", "S", "DS", "Hf", "SH", "Gf", "logKf"
        """

        match label:
            case "Cp":
                return self._data.Cp_list
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
