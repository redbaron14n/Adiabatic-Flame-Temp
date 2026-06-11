# ###################
# Ian Janes
# Prof. Don Lipkin
# Adiabatic Flame Temperature Calculator
# ###################

from domain.compound import Compound

class CombustionReaction:

    def __init__(
            self,
            fuels: dict[Compound, float],
            oxidants: dict[Compound, float],
            pressure: float,
            temps: dict[Compound, float]
        ):
        
        self._set_fuels(fuels)
        self._set_oxidants(oxidants)
        self.pressure = pressure
        self.temperatures = temps


    ########################################
    # Getters and Setters
    ########################################


    def _set_fuels(self, fuels: dict[Compound, float]):

        for amount in fuels.values():
            if amount <= 0:
                raise ValueError("Fuel amounts must be positive values.")
        self._fuels = self._normalize_ratio(fuels)


    def _set_oxidants(self, oxidants: dict[Compound, float]):

        for amount in oxidants.values():
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
    def temperatures(self, value: dict[Compound, float]):

        for temp in value.values():
            if temp < 0:
                raise ValueError("Temperature must be a non-negative value.")
        self._temperatures = value


    ########################################
    # Private Methods
    ########################################


    def _normalize_ratio(self, ratio: dict[Compound, float]) -> dict[Compound, float]:

        total = sum(ratio.values())
        return {compound: amount / total for compound, amount in ratio.items()}