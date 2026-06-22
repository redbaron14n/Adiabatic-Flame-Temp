from domain.dissociation import Dissociation
from domain.compounds import compounds
from numpy.typing import NDArray
from scipy.optimize import fsolve, least_squares, OptimizeResult
import numpy as np

INCLUDED: list[str] = sorted({
    "Methane",
    "Oxygen",
    "Carbon_Dioxide",
    "Water",
    "Carbon",
    "Carbon_Monoxide",
    "Hydrogen",
    "Hydrogen_Monatomic",
    "Hydroxyl",
    # "Oxygen_Monatomic"
})
print(f"INCLUDED: {INCLUDED}")
INIT_ATOMS: dict[int, float] = {1: 1.332, 6: 0.333, 8: 1.334}
INITIAL_ENTHALPY: float = -24.932709
MASS_BALANCE_WEIGHT: float = 1000.0
EQUILIBRIUM_WEIGHT: float = 1.0
MIN_LOG_GUESS: float = -50.0
ENERGY_WEIGHT: float = 1.0


def get_non_elementals() -> list[str]:

    non_elementals: list[str] = []
    for c_id in INCLUDED:
        if compounds[c_id].composition:
            non_elementals.append(c_id)
    # print(f"NON-ELEMENTALS: {non_elementals}")
    return non_elementals


def get_dissociation_objects(non_elementals: list[str]) -> dict[str, Dissociation]:

    dissociation_objects: dict[str, Dissociation] = {
        compound: Dissociation(compound, set(compounds[compound].composition.keys()))
        for compound in non_elementals
    }
    # print(f"DISSOCIATION_OBJECTS: {dissociation_objects}")
    return dissociation_objects


def get_index_dict() -> dict[str, int]:

    i = {c_id: i for i, c_id in enumerate(INCLUDED)}
    i["T"] = len(INCLUDED)
    # print(f"INDEX: {i}")
    return i


def get_present_atoms() -> list[int]:

    atoms: set[int] = set()
    for c_id in INCLUDED:
        atoms.update(compounds[c_id].atomic_composition().keys())
    present_atoms: list[int] = sorted(atoms)
    # print(f"PRESENT_ATOMS: {present_atoms}")
    return present_atoms


def get_residual_index_dict(present_atoms: list[int], non_elementals: list[str]) -> dict[str, int]:

    r: dict[str, int] = {}
    for indx, atom in enumerate(present_atoms):
        r[f"{atom}_mass_balance"] = indx
    offset = len(present_atoms)
    for c_id in non_elementals:
        r[c_id] = offset
        offset += 1
    r["Energy_balance"] = offset
    print(f"R: {r}")
    return r


def construct_init_log_guess(i: dict[str, int]) -> NDArray[np.float64]:

    init_log_guess: NDArray[np.float64] = np.full(len(INCLUDED) + 1, MIN_LOG_GUESS)
    init_log_guess[i["Carbon_Dioxide"]] = -0.47756
    init_log_guess[i["Oxygen"]] = -3
    init_log_guess[i["Water"]] = -0.17653
    init_log_guess[i["T"]] = 3000
    # print(f"INIT_LOG_GUESS: {init_log_guess}")
    return init_log_guess


def calc_mass_balance_residuals(log_guess: NDArray[np.float64], present_atoms: list[int], i: dict[str, int], r: dict[str, int]) -> NDArray[np.float64]:

    mbr = np.zeros_like(log_guess)
    for atom in present_atoms:
        guess_qty = 0
        for c_id in INCLUDED:
            guess_qty += compounds[c_id].atomic_composition().get(atom, 0) * (10**log_guess[i[c_id]])
        mbr[r[f"{atom}_mass_balance"]] = guess_qty - INIT_ATOMS[atom]
    return mbr


def calc_equil_residuals(log_guess: NDArray[np.float64], non_elementals: list[str], i: dict[str, int], r: dict[str, int], dissociation_objects: dict[str, Dissociation]) -> NDArray[np.float64]:

    er = np.zeros_like(log_guess)
    for c_id in non_elementals:
        dissociation = dissociation_objects[c_id]
        resid = dissociation.equilibrium_residual(log_guess, i)
        er[r[c_id]] = resid
    return er


def calc_energy_residual(log_guess: NDArray[np.float64], i: dict[str, int], r: dict[str, int]) -> float:

    product_enthalpy: float = 0
    temp = log_guess[i["T"]]
    for c_id, indx in i.items():
        if c_id == "T":
            continue
        compound = compounds[c_id]
        product_enthalpy += 10**log_guess[indx] * (compound.SH(temp) + compound.stdHf)
    energy_residual = product_enthalpy - INITIAL_ENTHALPY
    return energy_residual


def residual_function(log_guess: NDArray[np.float64], present_atoms: list[int], item_indices: dict[str, int], residual_indices: dict[str, int], non_elementals: list[str], dissociation_objects: dict[str, Dissociation]) -> NDArray[np.float64]:
    
    residuals = np.zeros_like(log_guess)
    residuals += MASS_BALANCE_WEIGHT * calc_mass_balance_residuals(log_guess, present_atoms, item_indices, residual_indices)
    residuals += EQUILIBRIUM_WEIGHT * calc_equil_residuals(log_guess, non_elementals, item_indices, residual_indices, dissociation_objects)
    residuals[residual_indices["Energy_balance"]] = ENERGY_WEIGHT * calc_energy_residual(log_guess, item_indices, residual_indices)
    return residuals


def main():

    index_dict = get_index_dict()
    init_log_guess = construct_init_log_guess(index_dict)
    present_atoms = get_present_atoms()
    non_elementals = get_non_elementals()
    residual_index_dict = get_residual_index_dict(present_atoms, non_elementals)
    dissociation_objects = get_dissociation_objects(non_elementals)
    solution = least_squares(
        residual_function,
        init_log_guess,
        args=(present_atoms, index_dict, residual_index_dict, non_elementals, dissociation_objects),
        max_nfev=10000,
        xtol=1e-12
    )
    return solution


print("Running main()...")
print(main())