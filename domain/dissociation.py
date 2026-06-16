# ###################
# Ian Janes
# Professor Don Lipkin
# Adiabatic Flame Temperature
# Dissociation Reaction Class File
# ###################

from chempy import balance_stoichiometry
from domain.compounds import compounds, compounds_by_formula
from math import isclose, log10
import numpy as np
from numpy.typing import NDArray

class Dissociation:

    def __init__(self, m_id: str, r_ids: set[str]):

        self._set_molecule(m_id)
        self._set_radicals(r_ids)
        self._set_stoichioemetry()
        self._set_nonsolids()


    ########################################
    # Getters and Setters
    ########################################


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
        rad_raw, prod_raw = balance_stoichiometry({mfrm}, rfrms)
        rad = {species: float(coeff) for species, coeff in rad_raw.items()} # Sanitizing result from SymPy objects
        prod = {species: float(coeff) for species, coeff in prod_raw.items()}
        stoich_dict = rad | prod
        factor = stoich_dict[mfrm]
        self._stoich: dict[str, float] = {k: v/factor for k, v in stoich_dict.items()}


    def _set_nonsolids(self):

        self._nonsolids: set[str] = set()
        for species in self._stoich.keys():
            compound = compounds_by_formula[species]
            if compound.state != "s":
                self._nonsolids.add(compound.id)


    ########################################
    # Private Methods
    ########################################


    def _validate_guess(self, log_guess: NDArray[np.float64], species_indices: dict[str, int]):

        if len(log_guess) != len(species_indices):
            raise ValueError(
                f"Guess list length does not match number of species plus 1.\n"
                f"Guess: {log_guess}\n"
                f"Species Indices: {species_indices}"
            )
        elif not all(species in species_indices for species in ({self._molecule} | self._radicals)):
            raise ValueError(
                f"Guess list does not contain all species in the reaction.\n"
                f"Guess: {log_guess}\n"
                f"Species Indices: {species_indices}\n"
                f"Reaction Species: {self._stoich.keys()}"
            )


    def _calc_pressure_exp(self) -> float:

        """
        Calculates and returns the pressure exponent for the equilibrium residual calculation based on the stoichiometry of the reaction.
        """

        stoich_dict = self._stoich
        mfrm = compounds[self._molecule].formula
        rfrms = {compounds[r].formula for r in (self._radicals & self._nonsolids)}
        return stoich_dict[mfrm] - sum(stoich_dict[r] for r in rfrms)
    

    def _calc_gas_moles(self, log_guess: NDArray[np.float64], species_indices: dict[str, int]) -> float:

        total_moles = 0.0
        for species, indx in species_indices.items():
            if (species != "T") and (compounds[species].state != "s"):
                total_moles += 10**log_guess[indx]
        return total_moles
        

    def _calc_log_conc_product(self, log_guess: NDArray[np.float64], species_indices: dict[str, int]) -> float:

        product = log_guess[species_indices[self._molecule]] * self._stoich[compounds[self._molecule].formula]
        for species, coeff in self._stoich.items():
            species = compounds_by_formula[species].id
            if (species in self._nonsolids) and (species != self._molecule):
                product -= log_guess[species_indices[species]] * coeff
        return product
    

    def _calc_log_pres_factor(self, log_guess: NDArray[np.float64], species_indices: dict[str, int], pressure: float) -> float:

        exponent = self._calc_pressure_exp()
        if isclose(exponent, 0.0):
            return 1.0
        fraction = pressure / self._calc_gas_moles(log_guess, species_indices)
        return exponent * log10(fraction)


    ########################################
    # Public Methods
    ########################################
    

    def get_log_eq_constant(self, temperature: float) -> float:

        compound = compounds[self._molecule]
        log_eq_constant = compound.logKf(temperature)
        return log_eq_constant


    def equilibrium_residual(self, log_guess: NDArray[np.float64], species_indices: dict[str, int], pressure: float = 1e5) -> float:

        self._validate_guess(log_guess, species_indices)
        temp = log_guess[species_indices["T"]]
        log_conc_product = self._calc_log_conc_product(log_guess, species_indices)
        log_pressure_factor = self._calc_log_pres_factor(log_guess, species_indices, pressure)
        log_ecc = self.get_log_eq_constant(temp)
        return log_conc_product + log_pressure_factor - log_ecc