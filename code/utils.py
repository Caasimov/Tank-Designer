import math as m
from typing import Tuple

def findt_1(P: float, R: float, sigma_y: float) -> float:
    '''
    Determine the value of the shell thickness t_1 due to hoop stress
    '''
    t_1 = P * R / (sigma_y * m.pow(10, 6))
    
    return t_1
    
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

def findGeoProperties(t_1: float, R: float) -> Tuple[float, float]:
    '''
    Determine the geometric properties of the structure
    '''
    I = m.pi * t_1 * R**3 # Determine area moment of inertia
    A = m.pi * ((R + (0.5 * t_1))**2 - (R - (0.5*t_1))**2) # Determine cross-sectional area
    
    return I, A

def eulerColBuck(F: float, E: float, l: float, I: float, A: float) -> Tuple[float, float, bool]:
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