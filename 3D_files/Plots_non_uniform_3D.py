import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
from scipy.special import erf, erfc
from scipy import optimize
from textwrap import wrap

plt.ioff()

NR=200
NP=200
NT=200

#Plotting optical depth for non-uniform grid
flat_t = []
den = [] * (NR+1 * NP+1 * NT+1)

dataset = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/text_files/density_NR=%s.txt' % NR,'r')

flat_t = np.loadtxt(dataset)

ThreeDt = np.reshape(flat_t,(NR+1,NP+1,NT+1),order='C')
print(np.shape(ThreeDt))

plt.figure()
plt.contourf(ThreeDt[:,:,-1])
plt.colorbar()
plt.title("Colourmap showing final optical depth as a function of theta and phi")
plt.xlabel("phi")
plt.ylabel("theta")
plt.show()
