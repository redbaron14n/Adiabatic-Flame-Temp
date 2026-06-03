# ###################
# Ian Janes
# Prof. Don Lipkin
# Adiabatic Flame Temperature
# Dissociative Flame Temperature Calculation Class File
# ###################

from compounds import compounds

class DissociativeReaction:

    def __init__(self, fuels: dict[str, float], oxi: dict[str, float], temps: dict[str, float], conc_res: int = 100, comb: float = 1.0):

        """
        :param dict[str, float] fuels: A dictionary mapping fuel species' IDs to their relative ratios.
        :param dict[str, float] oxi: A dictionary mapping oxidant species' IDs to their relative ratios.
        :param dict[str, float] temps: A mapping of species IDs to their corresponding temperatures in Kelvin for this reaction.
        :param int conc_res: The concentration resolution for the calculation.
        :param float comb: The percent combustion.
        """

        self._set_fuels(fuels)
        self._set_oxidants(oxi)
        self.temperatures = temps
        self.concentration_resolution = conc_res
        self.percent_combustion = comb


    @property
    def fuels(self) -> dict[str, float]:

        """
        :return: A dictionary mapping fuel species' IDs to their relative ratios.
        """

        return self._fuels
    

    def _set_fuels(self, fuels: dict[str, float]):

        if not all(species in compounds for species in fuels.keys()):
            raise ValueError("All species in fuels must be valid compound IDs.")
        if not all(ratio > 0 for ratio in fuels.values()):
            raise ValueError("All fuel ratios must be positive.")
        self._fuels: dict[str, float] = self._normalize_ratios(fuels)


    @property
    def oxidants(self) -> dict[str, float]:

        """
        :return: A dictionary mapping oxidant species' IDs to their relative ratios.
        """

        return self._oxidants
    

    def _set_oxidants(self, oxidants: dict[str, float]):

        if not all(species in compounds for species in oxidants.keys()):
            raise ValueError("All species in oxidants must be valid compound IDs.")
        if not all(ratio > 0 for ratio in oxidants.values()):
            raise ValueError("All oxidant ratios must be positive.")
        self._oxidants: dict[str, float] = self._normalize_ratios(oxidants)


    def _normalize_ratios(self, ratios: dict[str, float]) -> dict[str, float]:

        """
        Normalizes the given ratios so that they sum to 1 and returns the normalized ratios.

        :param dict[str, float] ratios: A dictionary mapping species' IDs to their relative ratios.

        :return: A dictionary mapping species' IDs to their normalized relative ratios.
        """

        total = sum(ratios.values())
        return {species: ratio / total for species, ratio in ratios.items()}


    @property
    def temperatures(self) -> dict[str, float]:

        """
        :return: A mapping of species IDs to their corresponding entry temperatures in Kelvin for this reaction.
        """

        return self._temps
    

    @temperatures.setter
    def temperatures(self, temps: dict[str, float]):

        if not all(species in self._fuels | self._oxidants for species in temps.keys()):
            raise ValueError("All species in temperatures must be in fuels or oxidants.")
        if not all(species in temps for species in self._fuels | self._oxidants):
            raise ValueError("All species in fuels and oxidants must have an entry in temperatures.")
        if not all(temp > 0 for temp in temps.values()):
            raise ValueError("All temperatures must be greater than 0 K.")
        self._temps = temps


    @property
    def concentration_resolution(self) -> int:

        """
        :return: How many ratios of fuels to oxidants to calculate for this reaction.
        """

        return self._conc_res
    

    @concentration_resolution.setter
    def concentration_resolution(self, res: int):

        if res < 2:
            raise ValueError("Concentration resolution must be 2 or greater.")
        self._conc_res = res


    @property
    def percent_combustion(self) -> float:

        """
        :return: The percent combustion to calculate the adiabatic flame temperature for.
        """

        return self._combustion
    

    @percent_combustion.setter
    def percent_combustion(self, comb: float):

        if not (0 < comb <= 1):
            raise ValueError("Percent combustion must be between 0 and 1.")
        self._combustion = comb