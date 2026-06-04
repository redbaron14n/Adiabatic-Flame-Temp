# ###################
# Ian Janes
# Prof. Don Lipkin
# Adiabatic Flame Temperature
# Dissociative Flame Temperature Calculation Class File
# ###################

from chempy.util.parsing import formula_to_composition
from domain.compounds import compounds
from numpy import column_stack, float64, linspace
from numpy.typing import NDArray

class DissociativeReaction:

    def __init__(self, fuels: dict[str, float], oxi: dict[str, float], temps: dict[str, float], conc_res: int = 100):

        """
        :param dict[str, float] fuels: A dictionary mapping fuel species' IDs to their relative ratios.
        :param dict[str, float] oxi: A dictionary mapping oxidant species' IDs to their relative ratios.
        :param dict[str, float] temps: A mapping of species IDs to their corresponding temperatures in Kelvin for this reaction.
        :param int conc_res: The concentration resolution for the calculation.
        """

        self._set_fuels(fuels)
        self._set_oxidants(oxi)
        self.temperatures = temps
        self.concentration_resolution = conc_res


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
        self._fuel_atoms = self._atoms_per_input(self._fuels)
        self._update_atom_counts()


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
        self._oxidant_atoms = self._atoms_per_input(self._oxidants)
        self._update_atom_counts()


    def _normalize_ratios(self, ratios: dict[str, float]) -> dict[str, float]:

        """
        Normalizes the given ratios so that they sum to 1 and returns the normalized ratios.

        :param dict[str, float] ratios: A dictionary mapping species' IDs to their relative ratios.

        :return: A dictionary mapping species' IDs to their normalized relative ratios.
        """

        total = sum(ratios.values())
        return {species: ratio / total for species, ratio in ratios.items()}
    

    def _atoms_per_input(self, ratio_dict: dict[str, float]) -> dict[int, float]:

        """
        :param dict[str, float] ratio_dict: A dictionary mapping species IDs to their relative ratios.

        :return: A dictionary mapping each atom present in the fuels and oxidants of this reaction to the total number of atoms contributed by the input ratios.
        """

        atom_totals: dict[int, float] = dict()
        for species in ratio_dict.keys():
            species_makeup = formula_to_composition(compounds[species].formula)
            for atom, count in species_makeup.items():
                atom_totals[atom] = atom_totals.get(atom, 0) + ratio_dict[species] * count
        return atom_totals


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
        self._update_atom_counts()


    def _update_atom_counts(self):

        """
        Updates the initial atoms counts for each fuel-to-oxidant ratio to reflect the percent combustion for this reaction.
        """

        required = (
            hasattr(self, "_fuel_atoms"),
            hasattr(self, "_oxidant_atoms"),
            hasattr(self, "_conc_res"),
        )

        if not all(required):
            return # Don't update if any of the required attributes haven't been set yet, as in initialization

        self._calc_fuel_oxi_ratios()
        self._calc_init_atom_counts()


    def _calc_fuel_oxi_ratios(self):

        """
        Generates the (n x 2) numpy array where each row is a combination of fuel and oxidant ratios to calculate for this reaction.
        """

        res = self._conc_res
        fuel_props = linspace(1/(res + 1), res/(res + 1), res)
        oxidant_props= 1.0 - fuel_props
        self._fuel_oxi_ratios = column_stack((fuel_props, oxidant_props))
    

    def _calc_init_atom_counts(self):

        """
        Generates the list of dictionaries mapping each atom present in the reactants to their total initial counts for each fuel-to-oxidant ratio to be calculated.
        """

        initial_atoms: list[dict[int, float64]] = []
        for ratio in self._fuel_oxi_ratios:
            atom_qtys: dict[int, float64] = dict()
            for atom, amount in self._fuel_atoms.items():
                atom_qtys[atom] = ratio[0] * amount # The atoms contributed by the fuels
            for atom, amount in self._oxidant_atoms.items():
                atom_qtys[atom] = atom_qtys.get(atom, 0) + ratio[1] * amount # Adding the atoms contributed by the oxidants
            initial_atoms.append(atom_qtys)
        self._init_atom_counts = initial_atoms