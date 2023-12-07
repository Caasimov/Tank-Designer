import json
import utils

# Load JSON file
with open('./data/materials.json', 'r') as fname:
    material_types = json.load(fname)

# Initialize variables
var_history = {}


# Recursion logic

for material_name, material_properties in material_types.items():

    var_history_buffer = [] # Reset the variable history buffer for the specific material under consideration
    num_iterations = 0
    
    while num_iterations <= 1000: #NOTE: REPLACE WITH ACTUAL CONDITION
        var_temp_buffer = [] # Reset the calculated design variable values
        
        
        var_history_buffer.append(var_temp_buffer)
        num_iterations += 1
    
    
    var_history[material_name] = var_history_buffer # Store iteration history in a material-specific dictionary

# Graphical output