from domain.reaction import Reaction
from numpy import array, vstack

test_reaction = Reaction({"Hydrogen", "Oxygen"}, {"Hydrogen": 298.15, "Oxygen": 298.15})

test_rows = [1, 12, 17, 23, 27, 28, 31, 34, 36, 38]
test_points = [x - 1 for x in test_rows]
test_concentrations = array([x / 51 for x in test_rows])

flame_table = test_reaction.calc_flame_table("Hydrogen", {"Hydrogen": 2, "Oxygen": 1}, 50)
test_values = array([flame_table[1, i] for i in test_points])

test_table = vstack((test_concentrations, test_values)).T
print(test_table)