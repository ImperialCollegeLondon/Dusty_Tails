import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt

T = np.linspace(1400, 4400, 1000)

alpha = 0.1
mu = 101.96
A = 77365
B = 39.3
density = 4.0

Mu = 1.67e-24
kB = 1.38e-16

s=10e-4

dt = []
for i in T:
    J = (-alpha/density)*np.exp(-(A/i)+B)*np.sqrt((mu*Mu)/(2*np.pi*kB))
    new = s/J
    dt.append(new)




plt.plot(T,dt)
plt.xlim(1000,5000)
plt.ylim(-7e8,1e8)
plt.show()
