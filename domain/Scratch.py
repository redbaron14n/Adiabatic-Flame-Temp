    def _calc_product_heat(self, guess: NDArray[float64]) -> float:

        temp: float = guess[self._item_indices["T"]]
        product_heat: float = 0.0
        for species, idx in self._item_indices.items():
            if species != "T":
                compound = compounds[species]
                product_heat += guess[idx] * (compound.SH(temp) + compound.stdHf)
        return product_heat
    

    def _calc_reactant_heat(self, fuel_ratio: float) -> float:

        for fuel, conc in self._fuels.items():
            temp = self.temperatures[fuel]
            compound = compounds[fuel]
            amount = conc * fuel_ratio


    def _energy_residual(self, guess: NDArray[float64]) -> float:

        product_heat = self._calc_product_heat(guess)


    def _residuals(self, guess: NDArray[float64], initial_atoms: dict[int, float64]) -> NDArray[float64]:

        """
        :param NDArray[float64] guess: The current guess for the species concentrations.
        :param dict[int, float64] initial_atoms: A dictionary mapping atom IDs to their initial counts.
        :param dict[str, int] species_indices: A mapping of species names to their corresponding indices in the guess list.

        :return: A numpy array of the residuals for the atom balances and dissociation equilibria for the given guess of species concentrations.
        """

        res_list: list[float] = []
        res_list.extend(self._abr_list(initial_atoms, guess))
        res_list.extend(self._er_list(guess))
        return array(res_list, dtype=float64)
    

    def calculate_equilibriums(self):

        equilibriums: list[NDArray[float64]] = []
        fuel_oxi_ratios = self._fuel_oxi_ratios
        guess_size = len(self._item_indices)
        for row in range(fuel_oxi_ratios.shape[0]):
            guess = ones(guess_size, dtype=float64)
            equil = fsolve(self._residuals, guess, args=(self._init_atom_counts[row]))
            print(equil)