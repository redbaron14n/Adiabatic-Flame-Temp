# ###################
# Ian Janes
# Professor Don Lipkin
# Adiabatic Flame Temperature
# Dissociation Reaction Class File
# ###################

from chempy import balance_stoichiometry
from domain.compound import Compound
from domain.compounds import compounds

class Dissociation:

    def __init__(self, molecule: str, radicals: set[str]):

        self._set_molecule(molecule)
        self._set_radicals(radicals)
        self._set_stoichioemetry()


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