# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# CSV Editor File
# ###################

import pandas as pd
from pathlib import Path
import numpy as np

MASTER_FILE = "thermochemical_data.csv"
REQUIRED_HEADER = ["Compound", "T", "Cf", "S", "(G-H)/T", "SH", "Hf", "G", "logKf"]

def get_csv_files() -> list[Path]:

    folder = Path()
    csv_files = []
    for file in folder.iterdir():
        if ((file.suffix == ".csv") and (file.name != "thermochemical_data.csv")):
            csv_files.append(file)
    return csv_files

def select_file(csv_files: list[Path]) -> Path | None:

    for i, file in enumerate(csv_files, 1):
        print(f"{i}: {file.name}")
    print("Q: Cancel")
    while True:
        choice = input("Select a file by number: ")
        if choice.lower() == "q":
            return None
        try:
            index = int(choice) - 1
            if 0 <= index < len(csv_files):
                return csv_files[index]
        except:
            pass
        print("Invalid input. Please enter a number or 'Q' to cancel.")

def validate_file(df: pd.DataFrame, file_name: str) -> bool:

    valid = True

    if list(df.columns) != REQUIRED_HEADER: # Checks header
        print(f"Invalid header in file {file_name}")
        valid = False
    
    if not df["Compound"].apply(lambda x: isinstance(x, str)).all(): # Checks if all entries in "Compound" column are strings
        bad_rows = df.index[~df["Compound"].apply(lambda x: isinstance(x, str))] # Finds all non-string entries
        for r in bad_rows:
            print(f"Invalid compound name in {file_name} at row {r+1}")
        valid = False
    
    numeric_cols = REQUIRED_HEADER[1:] # All but first column must be numeric
    for col in numeric_cols:
        for i, val in enumerate(df[col], start=2):
            if pd.isna(val):
                print(f"Invalid value (NaN) in {file_name} at row {i}, column {col}")
                valid = False
            elif not (np.isfinite(val) or np.isinf(val)):
                print(f"Invalid numeric value in {file_name} at row {i}, column {col}")
                valid = False

    return valid

def merge_csv():

    csv_files = get_csv_files()
    csv_file = select_file(csv_files)
    if csv_file is None:
        print("No file selected.")
        return
    
    new_df = pd.read_csv(csv_file)
    if not validate_file(new_df, csv_file.name):
        print("File validation failed.")
        return
    
    master = Path(MASTER_FILE)
    if master.exists():
        master_df = pd.read_csv(master)
    else:
        master_df = pd.DataFrame(columns=REQUIRED_HEADER) # Creates DataFrame with header to become master file later in process

    new_compounds = set(new_df["Compound"])
    conflicts = new_compounds.intersection(set(master_df["Compound"]))
    if conflicts:
        print(f"The following compounds already exist in the master file:")
        for c in conflicts:
            print(f" - {c}")
        confirm = input("Replace existing data? (y/n): ")
        if confirm != "y":
            print("Merge cancelled.")
            return
        master_df = master_df[~master_df["Compound"].isin(conflicts)]
    merged_df = pd.concat([master_df, new_df], ignore_index=True)
    merged_df.to_csv(master, index=False)
    print("Merge completed successfully.")

    delete = input(f"Delete {csv_file.name}? (y/n): ")
    if delete.lower() == "y":
        csv_file.unlink()
        print(f"{csv_file.name} deleted.")

merge_csv() # Temp call for testing