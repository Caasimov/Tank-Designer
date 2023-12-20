import math as m
import matplotlib.pyplot as plt
from typing import Tuple

def check_t_1(P: float, t_1: float, R: float, sigma_y: float) -> bool:
    '''
    Evaluate whether the pressure vessel withstands internal loads
    '''
    sigma = (P * R)/t_1
    
    return sigma < (sigma_y * m.pow(10, 6))
    
def findLambda(t_1: float, l: float, R: float, nu: float) -> float:
    '''
    Determine the value of constant lambda in the shell buckling equation
    '''
    base = R**2 * t_1**2 / (12 * (1 - nu**2))
    m_const = m.pi * l * m.pow(base, -0.25) # Determine the number of half-waves alond the length, l
    lambda_const = m_const * m.pi * R / l
    
    return lambda_const

def findk(t_1: float, l: float, R: float, nu: float) -> float:
    '''
    Determine the constant k in the shell buckling equation
    '''
    lambda_const = findLambda(t_1, l, R, nu)
    k = lambda_const + (12 / m.pi**4) * (l**4 / (R**2 * t_1**2)) * (1 - nu**2) / lambda_const
    
    return k

def findQ(P: float, t_1: float, E: float, R: float) -> float:
    '''
    Determine exponent Q in the shell buckling equation
    '''
    return (P / (E * m.pow(10, 9))) * (R / t_1)**2

def findl(V: float, R: float, t_1: float )-> float:
    '''Finding the required total length for a given tank radius and and the needed fuel volume'''
    return V/(m.pi*(R-t_1)**2) - 4*(R-t_1)/3

def roundUpN(val: float, N: int) -> float:
    '''Rounding up numbers to a required number of decimal places'''
    temp = val * 10**N
    return m.ceil(temp)/10**N

def findMass(R: float,l: float,t_1: float, rho: float)->float:
    '''find the mass of the empty tank from the given variables'''
    VS = 4*m.pi*(R**3-(R-t_1)**3)/3
    VC = m.pi*l*( R**2-(R-t_1)**2)
    return (VS + VC) * rho * m.pow(10, 3)

def findConnectorMass(R_outer: float, R_tank: float, F: float, E: float, rho: float, num_beam_pairs: int) -> Tuple[float, float]:
    '''
    Determine the connector properties of the tank
    '''
    l = m.sqrt(R_tank**2 + (R_outer - R_tank)**2) # Determine the lenngth of the hinge
    
    theta = m.atan(R_tank / (R_outer - R_tank))
    
    F_beam = F / (2 * num_beam_pairs * m.sin(theta))
    
    R = m.pow(4 * F_beam * l**2 / (m.pi**3 * E * m.pow(10, 9)), 0.25) # Determine connector radius
    
    mass = (l * m.pi * R**2) * rho * m.pow(10, 3)
    
    mass_total = mass * 2 * num_beam_pairs
    
    return mass_total, R
    
def findGeoProperties(t_1: float, R: float) -> Tuple[float, float]:
    '''
    Determine the geometric properties of the structure
    '''
    I = m.pi * t_1 * R**3 # Determine area moment of inertia
    A = m.pi * ((R)**2 - (R - t_1)**2) # Determine cross-sectional area

    
    return I, A

def eulerColBuckling(F: float, E: float, l: float, I: float, A: float) -> Tuple[float, float, bool]:
    '''
    Determine whether the structure will fail through Euler collumn buckling
    '''

    sigma_crit = (m.pi**2) * E * m.pow(10, 9) * I / (A*(l**2))

    sigma = F / A

    return sigma, sigma_crit, sigma < sigma_crit

def shellBuckling(P: float, F: float, k: float, t_1: float, l: float, R: float, A: float, E: float, nu: float) -> Tuple[float, float, bool]:
    '''
    Determine whether the structure will fail under shell buckling
    '''
    Q = findQ(P, t_1, E, R)
    
    sigma = F / A
    sigma_crit = (1.983 - (0.983 * m.exp(-23.14 * Q))) * k * (((m.pi**2) * E * m.pow(10, 9)) / (12 * (1 - nu**2))) * (t_1 / l)**2
    
    return sigma, sigma_crit, sigma < sigma_crit

def plotMaterialData(data, variable_properties):
    for material_name, connectors_data in data.items():
        # Determine the number of variables from the first entry
        num_variables = len(connectors_data[next(iter(connectors_data))][0])

        # Adjust here to exclude the last two variables
        num_variables_to_plot = num_variables - 3

        # Calculate the number of rows and columns for the subplots
        num_columns = 2
        num_rows = m.ceil(num_variables_to_plot / num_columns)

        # Adjust the figure size here if needed
        fig_width = 6
        fig_height = num_rows * 2  # Adjust the height as needed

        # Create a new figure for each material
        fig, axs = plt.subplots(num_rows, num_columns, figsize=(fig_width, fig_height))
        fig.canvas.manager.set_window_title(f'Material Analysis for {material_name}')
        axs = axs.flatten()  # Flatten the array to make indexing easier

        # Iterate over each variable, except the last two
        for var_index in range(num_variables_to_plot):
            ax = axs[var_index]
            ax.set_title(f'{variable_properties[0][var_index]}, {variable_properties[1][var_index]}')

            # Plot data for each number of connectors
            for num_connectors, histories in connectors_data.items():
                iterations = range(len(histories))
                ax.plot(iterations, [history[var_index] for history in histories], label=f'{num_connectors}')

            ax.legend()
            ax.set_xlabel('Iterations')
            ax.set_ylabel(f'{variable_properties[1][var_index]} {variable_properties[2][var_index]}')

        # Adjust layout to prevent overlap
        plt.tight_layout()

        plt.show()