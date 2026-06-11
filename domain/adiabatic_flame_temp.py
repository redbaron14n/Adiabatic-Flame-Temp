# ###################
# Ian Janes
# Prof. Don Lipkin
# Adiabatic Flame Temperature Calculator
# ###################

from domain.compounds import compounds

class CombustionReaction:

    def __init__(
            self,
            fuels: dict[str, float],
            oxidants: dict[str, float],
            temps: dict[str, float],
            pressure: float = 1e5,
            conc_res: int = 100
        ):

        """
        Initializes a new instance of the CombustionReaction class.

        :param dict[str, float] fuels: A dictionary mapping fuel compounds to their respective amounts.
        :param dict[str, float] oxidants: A dictionary mapping oxidant compounds to their respective amounts.
        :param dict[str, float] temps: A dictionary mapping compounds to their respective temperatures [K].
        :param float pressure: The pressure of the combustion reaction [Pa].
        :param int conc_res: The concentration resolution for the calculation.
        """
        
        self._set_fuels(fuels)
        self._set_oxidants(oxidants)
        self.pressure = pressure
        self.temperatures = temps
        self.concentration_resolution = conc_res
        self._set_conc_list()


    ########################################
    # Getters and Setters
    ########################################


    def _set_fuels(self, fuels: dict[str, float]):

        for fuel, amount in fuels.items():
            if fuel not in compounds:
                raise ValueError("Invalid fuel compound.")
            if amount <= 0:
                raise ValueError("Fuel amounts must be positive values.")
        self._fuels = self._normalize_ratio(fuels)


    def _set_oxidants(self, oxidants: dict[str, float]):

        for oxidant, amount in oxidants.items():
            if oxidant not in compounds:
                raise ValueError("Invalid oxidant compound.")
            if amount <= 0:
                raise ValueError("Oxidant amounts must be positive values.")
        self._oxidants = self._normalize_ratio(oxidants)

    
    @property
    def pressure(self):

        return self._pressure
    

    @pressure.setter
    def pressure(self, value: float):

        if value <= 0:
            raise ValueError("Pressure must be a positive value.")
        self._pressure = value


    @property
    def temperatures(self):

        return self._temperatures
    

    @temperatures.setter
    def temperatures(self, value: dict[str, float]):

        for temp in value.values():
            if temp < 0:
                raise ValueError("Temperature must be a non-negative value.")
        self._temperatures = value


    @property
    def concentration_resolution(self):

        return self._conc_res
    

    @concentration_resolution.setter
    def concentration_resolution(self, value: int):

        if value < 2:
            raise ValueError("Concentration resolution must be at least 2.")
        self._conc_res = value


    def _set_conc_list(self):

        denom = self._conc_res + 1 # Allows for desired resolution without 0 and 1 ratios (no flame)
        conc_list: list[dict[str, float]] = []
        for i in range(1, denom):
            fuel_ratio = i / denom
            conc_list.append({fuel: fuel_mix*fuel_ratio for fuel, fuel_mix in self._fuels.items()})
            conc_list[i-1].update({oxidant: oxi_mix*(1-fuel_ratio) for oxidant, oxi_mix in self._oxidants.items()})
        self._conc_list = conc_list


    ########################################
    # Private Methods
    ########################################


    def _normalize_ratio(self, ratio: dict[str, float]) -> dict[str, float]:

        total = sum(ratio.values())
        return {compound: amount / total for compound, amount in ratio.items()}
    

    def _calc_init_atoms(self, conc: dict[str, float]) -> dict[int, float]:

        """
        For given fuel and oxidant concentrations, calculates the initial number of each type of atom in the mixture.
        """

        init_atoms: dict[int, float] = {}
        for compound_id, amount in conc.items():
            compound = compounds[compound_id]
            for atom, count in compound.atomic_composition().items():
                init_atoms[atom] = init_atoms.get(atom, 0) + amount * count
        return init_atoms