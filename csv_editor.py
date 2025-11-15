# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# CSV Editor File
# ###################

import csv
from pathlib import Path

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
    invalid = True
    while invalid:
        choice = input("Select a file by number: ")
        try:
            index = int(choice) - 1
            if 0 <= index < len(csv_files):
                return csv_files[index]
            else:
                print("Invalid input. Please enter a number or 'Q' to cancel.")
        except:
            if choice.lower() == "q":
                return None
            else:
                print("Invalid input. Please enter a number or 'Q' to cancel.")

def row_validation(file) -> bool: # Figure out class of open file

    reader = csv.reader(file)
    for row_index, row, in enumerate(reader, start=1):
        # WIP

def validate_file(file_path: Path) -> bool:

    valid_file = True
    with open(file_path, "r") as file:
        if file.readline().strip().split(",") != REQUIRED_HEADER:
            valid_file = False
        elif row_validation(file) == False:
            valid_file = False
    return valid_file

def main():

    csv_files = get_csv_files()
    print(csv_files)
    csv_file = select_file(csv_files)
    print(csv_files)

main()