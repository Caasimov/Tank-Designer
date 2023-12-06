import math as m

def findt(P: float, R: float, sigma_y: float) -> float:
    '''
    Determine the value of the shell thickness t_1 due to hoop stress
    '''
    t_1 = P * R / (sigma_y * m.pow(10, 6))
    return t_1
    
def findLambda(P: float, l: float, R: float, sigma_y: float, nu: float) -> float:
    '''
    Determine the value of constant lambda in the shell buckling equation
    '''
    t_1 = findt(P, R, sigma_y) # Determine the current thickness of the shell
    base = R**2 * t_1**2 / (12 * (1 - nu**2))
    m_const = m.pi * l * m.pow(base, -0.25) # Determine the number of half-waves alond the length, l
    lambda_const = m_const * m.pi * R / l
    
    return m, lambda_const

def findk(P: float, l: float, R: float, sigma_y: float, nu: float) -> float:
    '''
    Determine the constant k in the shell buckling equation
    '''
    t_1 = findt(P, R, sigma_y)
    _, lambda_const = findLambda(P, l, R, sigma_y, nu)
    k = lambda_const + (12 / m.pi**4) * (l**4 / (R**2 * t_1**2)) * (1 - nu**2) / lambda_const
    
    return k
    

def findQ(P: float, E: float, R: float, sigma_y: float) -> float:
    t = findt(P,R,sigma_y) #input from pressure function 
    return( (P/E)*(R/t)**2)

def eulerColBuck(F: float, E: float, l: float, R: float, P: float, sigma_y: float) -> bool:
    '''
    here we are determinining if the structure will fail through collumn buckling
    '''
    t = findt(P,R,sigma_y) #input from pressure function 
    I = m.pi * t * R**3
    A = m.pi *((R+0.5*t)**2-(R-0.5*t)**2)
    
    sigmaCrit = (m.pi**2 )* E *I /(A*(l**2))

    sigma = F/A

    if sigma > sigmaCrit:
        return False 

    else:
        return True