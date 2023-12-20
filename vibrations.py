import numpy as np
from scipy.linalg import eigvals
import math
import matplotlib.pyplot as plt

#Spacecraft information
m1 = 40.31 #kg
r1 = 0.7875 #m
t1 = 0.0015 #m
l1 = 1.9598 #m
E1 = 72 * pow(10,9) #Pa

#Tank information
m2 = 24.57 #kg
r2 = 0.43663 #m
t2 = 0.00061 #m
l2 = 1.0014 #m
E2 = 210 * pow(10, 9) #Pa


#calculate area of cylinder cross-section
def area(r, t):
    return 2 * math.pi * r * t

#calculate spring stiffness
def k(L, r, t, E):
    return (2 * area(r, t) * E) / L

k1 = k(l1, r1, t1, E1)
k2 = k(l2, r2, t2, E2)

#calculate natural frequency
def frequency(m, L, r, t, E):
    A = area(r, t)
    return (1 / (2 * math.pi)) * math.sqrt((2 * A * E) / (m * L))

#mass matrix
M = np.array([[m1, 0],
              [0, m2]])

#stiffness matrix
K = np.array([[k1 + k2, -k2],
              [-k2, k2]])

#calculate eigenvalues
eigenvalues = eigvals(K, M)

#calculate nat frequencies
natural_frequencies = np.sqrt(eigenvalues)

print("Natural Frequencies:", natural_frequencies)

g = 9.81 #m/s^2
n = 8
m = 24.57 #kg
wn = frequency(m, l2, r2, t2, E2) #Hz
wf = 100 #Hz

#time list
t = np.linspace(0, 0.2, 10000)

#particular sol
xp = (8 * g / (m * (wn ** 2 - wf ** 2))) * np.sin(wf * t)

#homogeneous sol
xh = (8 * g / (m * (wn ** 2 - wf ** 2))) * (wf / wn) * np.sin(wn * t)

#final sol
x = xh+xp

## plot displacement
# plt.ylabel('displacement [m]')
# plt.xlabel('time [s]')
# plt.plot(t, x)
# plt.savefig('amplitutde.png')
plt.tight_layout()

#forced frequency list from 0 to 100
wf1 = np.arange(0, 100, 1)

#amplitutde formula
a = 8*g/(m*(wn**2 - wf1**2))

#plot amplitude
plt.ylabel('amplitutde [m]')
plt.xlabel('$\u03C9	_f$ [Hz]')
plt.plot(wf1, a)
plt.savefig('amplitutde.png')
plt.show()

