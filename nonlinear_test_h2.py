from domain.compound import Compound
from domain.compounds import compounds
from domain.dissociation import Dissociation
from numpy.typing import NDArray
from scipy.optimize import least_squares
import numpy as np


i: dict[str, int] = {
    "Hydrogen": 0,
    "Hydrogen_Monatomic": 1,
    "Hydroxyl": 2,
    "Oxygen": 3,
    "Oxygen_Monatomic": 4,
    "Water": 5,
    "T": 6
}

r: dict[str, int] = {
    "H_mass_balance": 0,
    "O_mass_balance": 1,
    "Hydrogen_Monatomic": 2,
    "Hydroxyl": 3,
    "Oxygen_Monatomic": 4,
    "Water": 5,
    "Energy_balance": 6
}

init_log_guess: NDArray[np.float64] = np.array([-3, -50, -50, -50, -50, -0.17653, 3000])

init_atoms: dict[int, float] = {1: 1.334, 8: 0.666}

initial_enthalpy: float = 0


def mass_balance_residuals(log_guess: NDArray[np.float64]) -> NDArray[np.float64]:

    mbr = np.zeros_like(log_guess)
    guess_H = 2*(10**log_guess[i["Hydrogen"]]) + 10**log_guess[i["Hydrogen_Monatomic"]] + 10**log_guess[i["Hydroxyl"]] + 2*(10**log_guess[i["Water"]])
    mbr[r["H_mass_balance"]] = init_atoms[1] - guess_H
    guess_O = 10**log_guess[i["Hydroxyl"]] + 2*(10**log_guess[i["Oxygen"]]) + 10**log_guess[i["Oxygen_Monatomic"]] + 10**log_guess[i["Water"]]
    mbr[r["O_mass_balance"]] = init_atoms[8] - guess_O
    return mbr # Consider putting this in log-space


def equil_residuals(log_guess: NDArray[np.float64]) -> NDArray[np.float64]:

    ebr = np.zeros_like(log_guess)
    for compound in ["Hydrogen_Monatomic", "Hydroxyl", "Oxygen_Monatomic", "Water"]:
        diss_obj = Dissociation(compound, set(compounds[compound].composition.keys()))
        resid = diss_obj.equilibrium_residual(log_guess, i)
        ebr[r[compound]] = resid
    return ebr


def energy_residual(log_guess: NDArray[np.float64]) -> float:

    product_enthalpy: float = 0
    temp = log_guess[i["T"]]
    for compound_id, indx in i.items():
        if compound_id == "T":
            continue
        compound = compounds[compound_id]
        product_enthalpy += 10**log_guess[indx] * (compound.SH(temp) + compound.stdHf)
    return product_enthalpy - initial_enthalpy


def residual_function(log_guess: NDArray[np.float64]) -> NDArray[np.float64]:

    residuals = np.zeros_like(log_guess)
    residuals += mass_balance_residuals(log_guess)
    residuals += equil_residuals(log_guess)
    residuals[r["Energy_balance"]] = energy_residual(log_guess)
    return residuals


def equilibrate(init_log_guess: NDArray[np.float64]):
    result = least_squares(residual_function, init_log_guess)
    return result


print(equilibrate(init_log_guess))