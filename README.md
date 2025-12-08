# Adiabatic-Flame-Temp

## UV package manager

Install UV on your PC.

```
pip install uv
```

In project directory run

```
uv sync
```

to create a populate .venv folder with the libraries used by the project.
To activate the virtual environment, run

```
source .venv/bin/activate
```

To deactivate the virtual environment, run

```
deactivate
```

To select the venv as the workspace venv, use ctl-shift-p 'python: select interpreter' and select the .venv/scripts/python.exe.

To run the app

```
uv run python .\app.py
```

Go to IP address printed in terminal.

At web app, at the top left there is dropdown labeled "Graph Mode" with two options; "reaction flame temperature" and "compound data".

### Reaction Flame Temperature

Under dropdown labeled "Select Reactants", add any and all reactants of your reaction, including inert reactants if any.

Under dropdown labeled "Selected Controlled Reactant", select which reactant's concentration will be used as an independent variable (x-axis) for the flame temperature.

If there are only two reactants, you may ignore "Ratios of Other Reactants", otherwise, enter the integer ratio of the listed reactants (i.e., for air, Nitrogen: 78, Oxygen: 21, Argon: 1)

Press "Update Graph" button once you have made desired selections.

### Compound Data

Under dropdown labeled "Select Compound", select the compound whose data you wish to view.

Under dropdown labeled "Select Variable", select the property you wish to view the data of.



## Adding Compounds and Reactions

### Merging Data

To add compounds to the data pool, add the .csv file containing the thermochemical data to the project directory. The .csv file should have the following header row otherwise only numbers or the string "inf" as entries:

'''
Compound,T,Cf,S,(G-H)/T,SH,Hf,G,logKf
'''

To merge the file into the master .csv file, run the csv_editor.py script

'''
uv run python -m .\csv_editor.py
'''

Follow the prompts to merge or extract data with/from the thermochemical_data.csv file. Final prompt asks to delete the given file to merge, just to stay tidy.

### IMPORTANT ###

Adding data of a compound already present in thermochemical_data.csv will replace all existing data with new data. This prevents duplicate points and inconsistent data. If you want to add data, keeping the pre-existing, extract the data using csv_editor.py, add your data to the created file, and then merge that file.

### Adding to Compound Dictionary

Go to .\domain\compounds.py and the following to the bottom:

'''
compounds["AAA"] = Compound(
    name = "BBB",
    formula = "CCC",
    id = "AAA",
    data = load_compound_data("AAA")
)
'''

where "AAA" is how the compound is named in the .csv file (i.e "Carbon_Dioxide"), "BBB" is how the compound should appear in dropdowns or texts (i.e. "Carbon Dioxide"), and "CCC" is the chemical formula of the compound (i.e. "CO2").

If an inert substance was added (i.e. Nitrogen gas or Argon), it must also be added to the set constant 'INERTS' found in config.py.

### Adding to Reaction Catalogue

Because guessing the products from the reactants is surprisingly challenging to implement, the current version of this program requires products to be hardcoded in. This may change in the future.

To add a new reaction, go config.py and under the function 'products_from_reactants', add the following *above* 'else:':

'''
elif (active_reactants == {"AAA", BBB"}) and not dissociation:
    products = {"CCC", "DDD"}
'''

where "AAA" and "BBB" (or more or less) are the active reactants, and "CCC" and "DDD" (or more or less) are the desired products. The 'dissociation' variable is future-proofing for when dissociation gets implemented.



# Planned updates

## Dissociation

Reactions at high temperatures have a myriad of additional products stemming from compounds breaking apart, which serve to significantly lower flame temperature.

## Reactants' Initial Temperatures

The initial temperatures of reactants plays a role in the sensible heat of reactants, which will influence the final flame temperature. For now, reactants are considered to enter the system at 298.15K.

## Exporting Data / Saving Graphs

It would be pragmatic if this program allowed the user to download calculated flame tables or printouts of graphs.

## Multiple Graphs

It may be very helpful to be able to have multiple graphs up at once to compare how changes to settings influence flame temperature.