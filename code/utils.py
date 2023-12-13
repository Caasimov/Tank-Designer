import math as m
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
    return (P / E) * (R / t_1)**2

def findl(V: float, R: float, t_1: float )-> float:
    '''Finding the required total length for a given tank radius and and the needed fuel volume'''
    return V/(m.pi*(R-t_1)**2) - 4*(R-t_1)/3

def roundUpN(val: float, N: int) -> float:
    '''rounding up numbers to a required number of decimal places'''
    temp = val * 10**N
    return m.ceil(temp)/10**N

def findMass(R: float,l: float,t_1: float, rho: float)->float:
    '''find the mass of the empty tank from the given variables'''
    VS = 4*m.pi*(R**3-(R-t_1)**3)/3
    VC = m.pi*l*( R**2-(R-t_1)**2)
    return (VS + VC) * rho * m.pow(10, 3)

def findGeoProperties(t_1: float, R: float) -> Tuple[float, float]:
    '''
    Determine the geometric properties of the structure
    '''
    I = m.pi * t_1 * R**3 # Determine area moment of inertia
    A = m.pi * ((R + (0.5 * t_1))**2 - (R - (0.5*t_1))**2) # Determine cross-sectional area
    
    return I, A

def eulerColBuckling(F: float, E: float, l: float, I: float, A: float) -> Tuple[float, float, bool]:
    '''
    Determine whether the structure will fail through Euler collumn buckling
    '''

    sigma_crit = (m.pi**2) * E * I / (A*(l**2))

    sigma = F / A

    return sigma, sigma_crit, sigma < sigma_crit

def shellBuckling(P: float, F: float, t_1: float, l: float, R: float, A: float, E: float, nu: float) -> Tuple[float, float, float, bool]:
    '''
    Determine whether the structure will fail under shell buckling
    '''
    Q = findQ(P, t_1, E, R)
    k = findk(t_1, l, R, nu)
    
    sigma = F / A
    sigma_crit = (1.983 - (0.983 * m.exp(-23.14 * Q))) * k * (((m.pi**2) * E) / (12 * (1 - nu**2))) * (t_1 / l)**2
    
    return k, sigma, sigma_crit, sigma < sigma_crit