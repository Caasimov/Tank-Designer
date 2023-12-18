import json
import utils


# Load JSON file
with open('./data/materials.json', 'r') as fname:
    material_types = json.load(fname)
    
# Spacecraft initial variables
mass_wet = 955.707
mass_tank = 59.11
vehicle_radius = 0.786
L_max = 2
vol_prop = 0.945
pressure = 24 * utils.m.pow(10, 5)
axial_acceleration = 9 * 9.80665

# Initialize variables
var_history = {}

num_squares_R = 1000
num_squares_t = 20

num_connectors_min = 3
num_connectors_max = 6

tol = utils.m.pow(10, -5)
delta_R = 0
delta_t = 0

# Recursion logic

for material_name, material_properties in material_types.items():

    var_history_buffer_main = {} # Reset the variable history buffer for the specific material under consideration
    num_iterations = 0
    R_0, R_max = 10*tol, 0.59
    t_0, t_max = 10*tol, 0.01
    
    for num_connector in range(num_connectors_min, num_connectors_max + 1):
        var_history_buffer = []
        num_iterations = 0
        R_0, R_max = 10*tol, 0.59
        t_0, t_max = 10*tol, 0.01
        while num_iterations <= 100: 
            # Logic check to verify whether sufficient iterations have been performed
            if len(var_history_buffer) >= 2:
                if utils.roundUpN(var_history_buffer[-2][5], 5) - utils.roundUpN(var_history_buffer[-1][5], 5) == 0.0:
                    break
                
            var_opt = [] # Reset the calculated design variable values
            
            delta_R = (R_max - R_0) / num_squares_R
            delta_t = (t_max - t_0) / num_squares_t
            
            # Iterate through grid to find optimal solution at current resolution
            for n in range(1, num_squares_R + 1):
                R_n = (n - 1) * delta_R + (delta_R / 2) + R_0
                for m in range(1, num_squares_t + 1):
                    t_m = (m - 1) * delta_t + (delta_t / 2) + t_0

                    l = utils.findl(vol_prop, R_n, t_m) # Determine the length of the cylinder
                    L = l + (2 * R_n) # Determine the total length of the fuel tank

                    # Ensure that l > R so that the shell buckling relations hold
                    if l > R_n and L < L_max:
                        mass = utils.findMass(R_n, l, t_m, material_properties['density']) # Determine mass of current fuel tank geometry
                        F = axial_acceleration * (mass_wet - mass_tank + m) # Determine new force due to tank geometry
                        mass_connector = utils.findConnectorMass(vehicle_radius, R_n, F, material_types['Steel 4130']['elasticity_modulus'], material_types['Steel 4130']['density'], num_connector)
                        F = axial_acceleration * (mass_wet - mass_tank + m + mass_connector) # Determine new force due to tank geometry
                        
                        k = utils.findk(t_m, l, R_n, material_properties['poisson_ratio']) # Determine the variable k to be optimized in shell buckling equation
                        I, A = utils.findGeoProperties(t_m, R_n) # Determine geometric properties of current cross-section
                        
                        # Check for hoop stress
                        check_hoop = utils.check_t_1(pressure, t_m, R_n, material_properties['tensile_strength'])
                        
                        # Check for Euler column buckling and shell buckling of the empty tank

                        sigma_eul, sigma_crit_eul, check_eul = utils.eulerColBuckling(F, material_properties['elasticity_modulus'], l, I, A)
                        sigma_shell, sigma_crit_shell, check_shell = utils.shellBuckling(pressure, F, k, t_m, l, R_n, A, material_properties['elasticity_modulus'], material_properties['poisson_ratio'])

                        # If both buckling checks pass and the mass of the current geometry is smaller than the previous optimal solution, then this is now recorded
                        if not var_opt and check_eul and check_shell and check_hoop:

                            new_var_opt = [l, R_n, t_m, k, mass_connector, mass, (sigma_eul, sigma_crit_eul), (sigma_shell, sigma_crit_shell)]
                            var_opt = new_var_opt
                            
                            # Bounds of grid search are updated to match with current optimal solution
                            R_0 = (n - 1) * delta_R + R_0
                            R_max = R_0 + delta_R
                            t_0 = (m - 1) * delta_t + t_0
                            t_max = t_0 + delta_t
                        
                        elif check_eul and check_shell and check_hoop and (mass + mass_connector < var_opt[4] + var_opt[5]):
                            new_var_opt = [l, R_n, t_m, k, mass_connector, mass, (sigma_eul, sigma_crit_eul), (sigma_shell, sigma_crit_shell)]
                            var_opt = new_var_opt
                            
                            # Bounds of grid search are updated to match with current optimal solution
                            R_0 = (n - 1) * delta_R + R_0
                            R_max = R_0 + delta_R
                            t_0 = (m - 1) * delta_t + t_0
                            t_max = t_0 + delta_t
                        
                    
            if var_opt:
                var_history_buffer.append(var_opt)
            num_iterations += 1
        
        index_name = f"N = {num_connector}"    
        var_history_buffer_main[index_name] = var_history_buffer # Store iteration history in a connector-number specific dictionary

    var_history[material_name] = var_history_buffer_main # Store iteration history in a material-specific dictionary

# Graphical output
print(var_history)

var_properties = [["Length of Cylinder", "Radius", "Thickness", "Shell Buckling Constant", "Connector Mass (Est.)", "Tank Mass"], 
            ["$l$", "$R$", "$t_1$", "$k$", "$m_{con}$", "$m$"], 
            ["[m]", "[m]", "[m]", "[-]", "[kg]", "[kg]"]]

utils.plotMaterialData(var_history, var_properties)