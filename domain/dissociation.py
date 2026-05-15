# ###################
# Ian Janes
# Professor Don Lipkin
# Adiabatic Flame Temperature
# Dissociation Reaction Class File
# ###################

from chempy import balance_stoichiometry
from domain.compound import Compound
from domain.compounds import compounds
from math import isclose

class Dissociation:

    def __init__(self, molecule: str, radicals: set[str]):

        self._set_molecule(molecule)
        self._set_radicals(radicals)
        self._set_stoichioemetry()
        self._set_nonsolids()


    def _set_molecule(self, molecule: str):

        self.molecule: str = molecule


    def _set_radicals(self, radicals: set[str]):

        self.radicals: set[str] = radicals


    def _set_stoichioemetry(self):

        radical_strs = {compounds[r].formula for r in self.radicals}
        molecule_str = compounds[self.molecule].formula
        rad, prod = balance_stoichiometry({molecule_str}, radical_strs)
        stoich_dict = rad | prod
        factor = stoich_dict[molecule_str]
        self._stoichiometry: dict[str, float] = {k: v/factor for k, v in stoich_dict.items()}


    def _set_nonsolids(self):

        self._nonsolids: set[str] = set()
        for species in self._stoichiometry.keys():
            for compound in compounds.values():
                if compound.formula == species and compound.state != "s":
                    self._nonsolids.add(compound.id)


    def _validate_guess(self, guess: list[float], species_indices: dict[str, int]):

        if len(guess) != len(species_indices) + 1:
            raise ValueError(
                f"Guess list length does not match number of species plus 1.\n"
                f"Guess: {guess}\n"
                f"Species Indices: {species_indices}"
            )
        elif not all(species in species_indices for species in self._stoichiometry):
            raise ValueError(
                f"Guess list does not contain all species in the reaction.\n"
                f"Guess: {guess}\n"
                f"Species Indices: {species_indices}\n"
                f"Reaction Species: {self._stoichiometry.keys()}"
            )


    def _calculate_pressure_exponent(self) -> float:

        """
        Calculates and returns the pressure exponent for the equilibrium residual calculation based on the stoichiometry of the reaction.
        """

        stoich_dict = self._stoichiometry
        molecule_str = compounds[self.molecule].formula
        radical_strs = {compounds[r].formula for r in self.radicals}
        return stoich_dict[molecule_str] - sum(stoich_dict[r] for r in radical_strs)
    

    def _calculate_pressure_factor(self, guess: list[float], species_indices: dict[str, int], pressure_bar: float) -> float:

        """
        Calculates and returns the pressure factor for the equilibrium residual calculation based on the current guess and pressure.

        @param guess (list[float]) : The current guess for the species concentrations, where the last element is the temperature.
        @param species_indices (dict[str, int]) : A mapping of species names to their corresponding indices in the guess list.
        @param pressure_bar (float) : The total pressure in bars.
        """

        exponent = self._calculate_pressure_exponent()
        if isclose(exponent, 0.0):
            return 1.0
        fraction = pressure_bar / sum(guess[species_indices[s]] for s in self._nonsolids)
        return fraction ** exponent


    def equilibrium_residual(self, guess: list[float], species_indices: dict[str, int], pressure_bar: float = 1.0) -> float:

        self._validate_guess(guess, species_indices)
        pass