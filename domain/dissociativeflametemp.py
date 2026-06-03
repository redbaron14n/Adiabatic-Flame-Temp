# ###################
# Ian Janes
# Prof. Don Lipkin
# Adiabatic Flame Temperature
# Dissociative Flame Temperature Calculation Class File
# ###################

from chempy.util.parsing import formula_to_composition
from domain.compounds import compounds

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
        self._set_init_conc_dicts()
        self._set_total_atoms()


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


    def _list_atoms(self) -> set[int]:

        """
        :return: A set of all atoms present in the fuels and oxidants of this reaction.
        """

        atoms: set[int] = set()
        for species in self._fuels | self._oxidants:
            formula = compounds[species].formula
            composition: dict[int, int] = formula_to_composition(formula)
            for element in composition.keys():
                atoms.add(element)
        return atoms


    def _find_atom_contributions(self) -> dict[int, dict[str, int]]:

        """
        :return: A dictionary mapping each atom present in the fuels and oxidants of this reaction to a dictionary mapping each species containing that atom to the number of atoms contributed by that species.

            For example, if the reaction is CH4 + 2 O2, the return value would be:
            {6: {'Methane': 1}, 1: {'Methane': 4}, 8: {'Oxygen': 2}}
        """

        atom_contrib: dict[int, dict[str, int]] = dict()
        atom_set = self._list_atoms()
        for atom in atom_set:
            atom_contrib[atom] = dict()
            for species in self._fuels | self._oxidants:
                formula = compounds[species].formula
                composition: dict[int, int] = formula_to_composition(formula)
                if atom in composition:
                    atom_contrib[atom][species] = composition[atom]
        return atom_contrib


    def _set_init_conc_dicts(self):

        """
        Generates a list of dictionaries mapping species IDs to their initial concentrations for each fuel-to-oxidant ratio to be calculated.
        """

        conc_step = 1 / (self._conc_res + 1)
        conc_list: list[dict[str, float]] = []
        fuel_conc = conc_step
        while fuel_conc < 1:
            oxidant_conc = 1 - fuel_conc
            conc_dict = {species: fuel_conc * ratio for species, ratio in self._fuels.items()}
            for species, ratio in self._oxidants.items():
                conc_dict[species] = oxidant_conc * ratio
            conc_list.append(conc_dict)
            fuel_conc += conc_step
        self._init_conc_list = conc_list


    def _set_total_atoms(self):

        atom_contrib = self._find_atom_contributions()
        atom_totals_list: list[dict[int, float]] = []
        for conc_dict in self._init_conc_list:
            atom_totals: dict[int, float] = dict()
            for atom, species_contrib in atom_contrib.items():
                for species, contrib in species_contrib.items():
                    atom_totals[atom] = atom_totals.get(atom, 0) + conc_dict[species] * contrib
            atom_totals_list.append(atom_totals)
        self._atom_totals_list = atom_totals_list

test = DissociativeReaction({'Methane': 1, "Carbon_Dioxide": 2}, {"Oxygen": 21, "Nitrogen": 79}, {'Methane': 300, 'Oxygen': 300, "Nitrogen": 300, "Carbon_Dioxide": 300}, conc_res=10, comb=0.8)
print(test._atom_totals_list)