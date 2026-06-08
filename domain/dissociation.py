# ###################
# Ian Janes
# Professor Don Lipkin
# Adiabatic Flame Temperature
# Dissociation Reaction Class File
# ###################

from chempy import balance_stoichiometry
from domain.compound import Compound
from domain.compounds import compounds, compounds_by_formula
from math import isclose

class Dissociation:

    def __init__(self, m_id: str, r_ids: set[str]):

        self._set_molecule(m_id)
        self._set_radicals(r_ids)
        self._set_stoichioemetry()
        self._set_nonsolids()


    @property
    def molecule_id(self) -> str:

        """
        :return: The ID of the molecule that is dissociating in this reaction.
        """

        return self._molecule


    def _set_molecule(self, id: str):

        self._molecule: str = id


    @property
    def radicals_ids(self) -> set[str]:

        """
        :return: A set of IDs for the radicals produced in this dissociation reaction.
        """

        return self._radicals


    def _set_radicals(self, ids: set[str]):

        self._radicals: set[str] = ids


    def _set_stoichioemetry(self):

        rfrms = {compounds[r].formula for r in self._radicals}
        mfrm = compounds[self._molecule].formula
        rad, prod = balance_stoichiometry({mfrm}, rfrms)
        stoich_dict = rad | prod
        factor = stoich_dict[mfrm]
        self._stoich: dict[str, float] = {k: v/factor for k, v in stoich_dict.items()}


    def _set_nonsolids(self):

        self._nonsolids: set[str] = set()
        for species in self._stoich.keys():
            compound = compounds_by_formula[species]
            if compound.state != "s":
                self._nonsolids.add(compound.id)


    def _validate_guess(self, guess: list[float], species_indices: dict[str, int]):

        if len(guess) != len(species_indices):
            raise ValueError(
                f"Guess list length does not match number of species plus 1.\n"
                f"Guess: {guess}\n"
                f"Species Indices: {species_indices}"
            )
        elif not all(species in species_indices for species in ({self._molecule} | self._radicals)):
            raise ValueError(
                f"Guess list does not contain all species in the reaction.\n"
                f"Guess: {guess}\n"
                f"Species Indices: {species_indices}\n"
                f"Reaction Species: {self._stoich.keys()}"
            )


    def _calculate_pressure_exponent(self) -> float:

        """
        Calculates and returns the pressure exponent for the equilibrium residual calculation based on the stoichiometry of the reaction.
        """

        stoich_dict = self._stoich
        mfrm = compounds[self._molecule].formula
        rfrms = {compounds[r].formula for r in (self._radicals & self._nonsolids)}
        return stoich_dict[mfrm] - sum(stoich_dict[r] for r in rfrms)
    

    def _calculate_pressure_factor(self, guess: list[float], species_indices: dict[str, int], pressure_bar: float) -> float:

        """
        Calculates and returns the pressure factor for the equilibrium residual calculation based on the current guess and pressure.

        :param list[float] guess: The current guess for the species concentrations, where the last element is the temperature.
        :param dict[str, int] species_indices: A mapping of species names to their corresponding indices in the guess list.
        :param float pressure_bar: The total pressure in bars.
        """

        exponent = self._calculate_pressure_exponent()
        if isclose(exponent, 0.0):
            return 1.0
        fraction = pressure_bar / sum(guess[species_indices[s]] for s in self._nonsolids)
        return fraction ** exponent
    

    def _calculate_concentration_product(self, guess: list[float], species_indices: dict[str, int]) -> float:

        """
        Calculates and returns the concentration product for the equilibrium residual calculation based on the current guess.

        :param list[float] guess: The current guess for the species concentrations, where the last element is the temperature.
        :param dict[str, int] species_indices: A mapping of species names to their corresponding indices in the guess list.
        """

        product = guess[species_indices[self._molecule]] ** self._stoich[compounds[self._molecule].formula]
        for species, coeff in self._stoich.items():
            species = compounds_by_formula[species].id # Convert from formula to ID for indexing guess list
            if (species in self._nonsolids) and (species != self._molecule):
                product /= guess[species_indices[species]] ** coeff
        return product


    def _validate_temperature(self, temperature: float):

        if temperature < 0:
            raise ValueError(f"Temperature must be greater than or equal to 0 K. Given: {temperature} K")


    def get_equilibrium_constant(self, temperature: float) -> float:

        """
        Calculates and returns the equilibrium constant for the reaction at the given temperature.

        :param float temperature: The temperature in Kelvin.
        """

        molecule_compound = compounds[self._molecule]
        ecc = 10 ** molecule_compound.logKf(temperature)
        return ecc


    def equilibrium_residual(self, guess: list[float], species_indices: dict[str, int], pressure_bar: float = 1.0) -> float:

        """
        Calculates and returns the equilibrium residual for the current guess of species concentrations and temperature at the given pressure.

        :param list[float] guess: The current guess for the species concentrations and temperature.
        :param dict[str, int] species_indices: A mapping of species names to their corresponding indices in the guess list.
        :param float pressure_bar: The total pressure in bars (default is 1.0 bar).
        """

        self._validate_guess(guess, species_indices)
        temp = guess[species_indices["T"]]
        self._validate_temperature(temp)
        conc_prod = self._calculate_concentration_product(guess, species_indices)
        pressure_factor = self._calculate_pressure_factor(guess, species_indices, pressure_bar)
        ecc = self.get_equilibrium_constant(temp)
        return conc_prod * pressure_factor - ecc