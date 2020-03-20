import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
from scipy.special import erf, erfc
from scipy import optimize
from textwrap import wrap

plt.ioff()

NR=20
NP=20
NT=20

error = []

#Plotting optical depth for non-uniform grid
#dens = [] * (NR+1 * NP+1 * NT+1)

opticaldepth_ana = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/opticaldepth_NR=%s.txt' % NR,'r')
phi = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/phi.txt', 'r')
theta = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/theta.txt', 'r')
density_ana = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/density_NR=%s.txt' % NR, 'r')

flat_phi = np.loadtxt(phi)
flat_theta = np.loadtxt(theta)

t_ana = np.loadtxt(opticaldepth_ana)
d_ana = np.loadtxt(density_ana)
# print (np.size(d_num))
# print(np.size(d_ana))
# print (np.size(t_num))
# print(np.size(t_ana))

t_ana = t_ana / max(t_ana)
d_ana = d_ana / max(d_ana)

t_ana = np.reshape(t_ana,(NR+1,NP+1,NT+1),order='C')
d_ana = np.reshape(d_ana,(NR+1,NP+1,NT+1),order='C')

plt.figure(1)
plt.contourf(flat_phi, flat_theta, t_ana[:,:,-1])
plt.colorbar()
plt.title("Analytical optical depth as a function of theta and phi")
plt.xlabel("phi")
plt.ylabel("theta")
plt.show()

plt.figure(2)
plt.contourf(flat_phi, flat_theta, d_ana[:,:,-1])
plt.colorbar()
plt.title("Analytical density as a function of theta and phi")
plt.xlabel("phi")
plt.ylabel("theta")
plt.show()
