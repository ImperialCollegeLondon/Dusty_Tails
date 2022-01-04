import numpy as np
from scipy.special import wofz as faddeeva
import matplotlib.pyplot as plt

clight = 2.998e10

def oscillator(omega0, omegap, sigma, gamma, omega):
    a = complex(omega/np.sqrt(2) * np.sqrt( np.sqrt(1+(gamma/omega)**2)+1 ),
                omega/np.sqrt(2) * np.sqrt( np.sqrt(1+(gamma/omega)**2)-1 ))
    return np.sqrt(3.14)*complex(0,1)*omegap**2/(np.sqrt(8)*sigma*a) * (faddeeva((a-omega0)/(np.sqrt(2)*sigma)) + faddeeva((a+omega0)/(np.sqrt(2)*sigma)))

def epsilon_from_oscillators(omega0_arr, omegap_arr, sigma_arr, gamma, epsinf, omega):
    eps = epsinf
    for i in range(0,len(omega0_arr)):
        eps += oscillator(omega0_arr[i], omegap_arr[i], sigma_arr[i], gamma, omega)
    return eps
        
# Fit between omega = 100 and 1400, in steps of 4. All in cm^-1.
omega0_arr = [1100, 983, 715, 384]
omegap_arr = [329, 709, 305, 469]
sigma_arr = [47, 57, 76, 113]
gamma = 4
epsinf = 3.61

omega_arr = []
re_epsilon_arr = []
im_epsilon_arr = []
for omega in range(100, 1400, 1):
    omega_arr.append(omega)
    eps = epsilon_from_oscillators(omega0_arr, omegap_arr, sigma_arr, gamma, epsinf, omega)
    re_epsilon_arr.append(eps.real)
    im_epsilon_arr.append(eps.imag)
    
plt.plot(omega_arr, re_epsilon_arr, c='k')
plt.plot(omega_arr, im_epsilon_arr, c='r', ls='--')
plt.ylim(0,7.5)
plt.xlim(100,1400)

plt.show()