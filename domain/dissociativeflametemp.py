# ###################
# Ian Janes
# Prof. Don Lipkin
# Adiabatic Flame Temperature
# Dissociative Flame Temperature Calculation Class File
# ###################

from chempy.util.parsing import formula_to_composition
from domain.compounds import compounds
from domain.dissociation import Dissociation
from config import determine_basic_products, INERTS
from numpy import append, array, column_stack, float64, linspace, ones
from numpy.typing import NDArray
from scipy.optimize import least_squares, OptimizeResult

MASS_WEIGHT = 100

class DissociativeReaction:

    def __init__(
            self,
            fuels: dict[str, float],
            oxi: dict[str, float],
            temps: dict[str, float],
            pres_bar: float = 1.0,
            conc_res: int = 100
        ):

        """
        :param dict[str, float] fuels: A dictionary mapping fuel species' IDs to their relative ratios.
        :param dict[str, float] oxi: A dictionary mapping oxidant species' IDs to their relative ratios.
        :param dict[str, float] temps: A mapping of species IDs to their corresponding temperatures in Kelvin for this reaction.
        :param float pres_bar: The pressure in bars for this reaction. Default is 1.0 bar.
        :param int conc_res: The concentration resolution for the calculation.
        """

        self._set_fuels(fuels)
        self._set_oxidants(oxi)
        self._set_reactants()
        self._set_products()
        self._set_item_indices()
        self.temperatures = temps
        self.pressure_bar = pres_bar
        self.concentration_resolution = conc_res
        self._update_fuel_oxi_ratios()
        self._set_init_ratios()
        self._set_bounds()
        self._set_dissociations()
        self._set_atomic_comps()
        self._reset_residual_labels()
        self._record = []
        self._pass_num = 0


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
            # species_makeup = formula_to_composition(compounds[species].formula)
            species_makeup = { # Sanitizing result from SymPy objects
                int(atom): float(count)
                for atom, count in formula_to_composition(
                    compounds[species].formula
                ).items()
            }
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
    def pressure_bar(self) -> float:

        """
        :return: The pressure in bars for this reaction.
        """

        return self._pres_bar
    

    @pressure_bar.setter
    def pressure_bar(self, pres: float):

        if pres <= 0:
            raise ValueError("Pressure must be greater than 0 bar.")
        self._pres_bar = pres


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

        self._update_fuel_oxi_ratios()
        self._calc_init_atom_counts()


    def _update_fuel_oxi_ratios(self):

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


    def _set_init_ratios(self): # Very similar to _calc_init_atom_counts; probably a way to merge these into one

        """
        Calculates and sets the initial ratios of each species for each fuel-to-oxidant ratio to be calculated.
        """

        init_ratios: list[dict[str, float]] = []
        for ratio in self._fuel_oxi_ratios:
            ratio_dict: dict[str, float] = dict()
            for fuel, fuel_ratio in self._fuels.items():
                ratio_dict[fuel] = fuel_ratio * ratio[0]
            for oxi, oxi_ratio in self._oxidants.items():
                ratio_dict[oxi] = oxi_ratio * ratio[1]
            init_ratios.append(ratio_dict)
        self._init_ratios = init_ratios


    def _set_reactants(self):

        """
        Sets the reactants for this reaction as the union of the fuels and oxidants.
        """

        fuels = set(self.fuels.keys())
        oxidants = set(self.oxidants.keys())
        self._reactants: set[str] = fuels.union(oxidants)


    def _set_products(self):
        
        """
        Finds and sets the product, including any dissociation products, and reactant (in-case of incomplete combustion) species IDs for the reaction's set of reactants.
        """

        active_reactants = self._reactants - INERTS
        fundamentals = determine_basic_products(active_reactants).union(self._reactants)
        products = set(fundamentals)
        for species in fundamentals:
            products.update(compounds[species].dissociates)
        self._products = products


    def _set_item_indices(self):

        """
        Sets the mapping of temperature and species IDs to their corresponding indices in the guess list for this reaction.
        """

        species_list = sorted(self._products) # Sort to ensure consistent ordering
        self._item_indices = {species: idx for idx, species in enumerate(species_list)}
        self._item_indices["T"] = len(species_list) # Temperature is the last element in the guess list


    def _set_atomic_comps(self):

        self._atomic_comps: dict[str, dict[int, float]] = dict()
        for species in self._products:
            self._atomic_comps[species] = {
                int(atom): float(count)
                for atom, count in formula_to_composition(
                    compounds[species].formula
                ).items()
            }


    def _atom_balance_residual(self, atom: int, initial_count: float, guess: NDArray[float64]) -> float:

        """
        :return: The residual for the balance of the given atom for the given guess of species concentrations, calculated as the initial count of the atom minus the count of the atom in the products based on the guess.
        """
        
        residual: float = -initial_count
        for species in self._products:
            atomic_comp = self._atomic_comps[species]
            residual += atomic_comp.get(atom, 0) * 10**(guess[self._item_indices[species]])
        return residual
    

    def _reset_residual_labels(self):

        self._residual_labels: list[str] = []


    def _abr_list(self, initial_atoms: dict[int, float64], guess: NDArray[float64]) -> list[float]:

        """
        :return: A list of the atom balance residuals for each atom present in the reactants for the given guess of species concentrations.
        """

        residual_list: list[float] = []
        for atom, count in initial_atoms.items():
            residual = self._atom_balance_residual(atom, count, guess)
            residual_list.append(residual)
            self._residual_labels.append(f"Atom {atom} balance")
        return residual_list
    

    def _set_dissociations(self):

        self._dissociations: dict[str, Dissociation] = dict()
        for molecule in self._products:
            if compounds[molecule].composition: # Only set dissociation reactions for non-elemental species
                components = set(compounds[molecule].composition.keys())
                self._dissociations[molecule] = Dissociation(molecule, components)
    

    def _eqr_list(self, guess: NDArray[float64]) -> list[float]:

        residuals: list[float] = []
        for diss_reaction in self._dissociations.values():
            residual = diss_reaction.equilibrium_residual(guess.tolist(), self._item_indices, self._pres_bar)
            residuals.append(residual)
            self._residual_labels.append(f"{diss_reaction._molecule} eq")
        return residuals
    

    def _calc_product_heat(self, guess: NDArray[float64]) -> float:

        temp: float = guess[self._item_indices["T"]]
        product_heat: float = 0.0
        for species, idx in self._item_indices.items():
            if species != "T":
                compound = compounds[species]
                product_heat += 10**guess[idx] * (compound.SH(temp) + compound.stdHf)
        return product_heat
    

    def _calc_reactant_heat(self, init_ratios: dict[str, float]) -> float:

        reactant_heat: float = 0.0
        for species, ratio in init_ratios.items():
            temp = self.temperatures[species]
            compound = compounds[species]
            reactant_heat += ratio * (compound.SH(temp) + compound.stdHf)
        return reactant_heat
    

    def _energy_residual(self, guess: NDArray[float64], init_ratios: dict[str, float]) -> float:

        """
        :param  NDArray[float64] guess: The current guess for the species concentrations and temperature.
        :param dict[str, float] init_ratios: A dictionary mapping species IDs to their initial ratios for this fuel-to-oxidant ratio.

        :return: The residual for the energy balance of the reaction for the given guess of species concentrations and temperature, calculated as the total enthalpy of the products based on the guess concentrations and temperature minus the total enthalpy of the reactants based on the initial ratios and entry temperatures.
        """

        product_heat = self._calc_product_heat(guess)
        reactant_heat = self._calc_reactant_heat(init_ratios)
        return product_heat - reactant_heat
    

    def _residual_list(self, guess: NDArray[float64], conc_indx: int) -> list[float]:

        self._reset_residual_labels()
        res_list: list[float] = []
        res_list.extend(self._abr_list(self._init_atom_counts[conc_indx], guess))
        res_list.extend(self._eqr_list(guess))
        res_list.append(self._energy_residual(guess, self._init_ratios[conc_indx]))
        self._residual_labels.append("Energy balance")

        if self._pass_num % 100 == 0:
            print(f"Pass {self._pass_num}")
            print(f"Guess: {guess}")
            print(res_list)
        self._pass_num += 1

        return res_list
    

    def _set_bounds(self):

        num_items = len(self._item_indices)
        self._lower_bounds = [-99] * (num_items - 1) + [0.0]
        self._upper_bounds = [5] * (num_items - 1) + [6000.0]
    

    def equilibrate(self, conc_indx: int) -> OptimizeResult:

        init_guess = [-0.3] * (len(self._item_indices) - 1) + [3000.0] # Initial guess of 0.5 for all species concentrations and 3000 K for temperature; can be adjusted as needed
        equil: OptimizeResult = least_squares(self._residual_list, init_guess, bounds=(self._lower_bounds, self._upper_bounds), args=(conc_indx,), ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=5000)
        
        for species, ratio in self._init_ratios[conc_indx].items():
            print(f"{species} initial amount: {ratio}")
        for species, idx in self._item_indices.items():
            if species != "T":
                print(f"{species} final amount: {10**equil.x[idx]}")
        print(f"Final temperature: {equil.x[self._item_indices['T']]} K")

        for label, residual in zip(self._residual_labels, equil.fun):
            print(f"{label} residual: {residual}")

        return equil