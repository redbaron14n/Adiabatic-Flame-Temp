# ###################
# Ian Janes
# Professor Don Lipkin
# MSEN 210 200
# Adiabatic Flame Temperature
# CSV Editor File
# ###################

from pathlib import Path

FOLDER_PATH = "newdata"
REQUIRED_HEADER = ["Compound", "T", "Cf", "S", "(G-H)/T", "SH", "Hf", "G", "logKf"]

def get_csv_files(path: str) -> list[Path]:

    folder = Path(path)
    csv_files = [f for f  in folder.iterdir() if f.suffix == ".csv"]
    return csv_files

def filter_csv_files(files: list[Path], required_header: list[str]) -> list[Path]:

    filtered_files = []
    for file in files:
        with open(file, 'r') as f:
            header = f.readline().strip().split(',')
            if header == required_header:
                filtered_files.append(file)
            else:
                print(f"File {file} skipped due to incorrect header.")
    return filtered_files

def main():

    csv_files = get_csv_files(FOLDER_PATH)
    print(csv_files)
    csv_files = filter_csv_files(csv_files, REQUIRED_HEADER)
    print(csv_files)

main()