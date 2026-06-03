# ###################
# Ian Janes
# Prof. Don Lipkin
# Adiabatic Flame Temperature
# Dissociative Flame Temperature Calculation Class File
# ###################

class DissociativeReaction:

    def __init__(self, fuels: set[str], oxi: set[str], temps: dict[str, float], conc_res: int = 100):

        """
        :param set[str] fuels: A set of fuel species IDs involved in this dissociative reaction.
        :param set[str] oxi: A set of oxidant species IDs involved in this dissociative reaction.
        :param dict[str, float] temps: A mapping of species IDs to their corresponding temperatures in Kelvin for this reaction.
        :param int conc_res: The concentration resolution for the calculation.
        """

        self._set_fuels(fuels)
        self._set_oxidants(oxi)
        self._temps = temps
        self._conc_res = conc_res


    @property
    def fuels(self) -> set[str]:

        """
        :return: A set of fuel species IDs involved in this dissociative reaction.
        """

        return self._fuels
    

    def _set_fuels(self, fuels: set[str]):

        self._fuels: set[str] = fuels


    @property
    def oxidants(self) -> set[str]:

        """
        :return: A set of oxidant species IDs involved in this dissociative reaction.
        """

        return self._oxidants
    

    def _set_oxidants(self, oxidants: set[str]):

        self._oxidants: set[str] = oxidants


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