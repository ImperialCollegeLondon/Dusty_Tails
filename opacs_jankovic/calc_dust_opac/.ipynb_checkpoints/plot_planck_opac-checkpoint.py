import numpy as np
import matplotlib.pyplot as plt



directory = 'corundum_K95'
#directory = 'enstatite_J98_J94_D95'
#directory = 'fayalite_F01'
#directory = 'graphite_D84'
#directory = 'olivine_F01'
#directory = 'silicon_carbide_L93'

sizes = np.loadtxt(directory+'/opac_sizes.dat')
temp = np.loadtxt(directory+'/opac_temp.dat')
label = '_0.1'
planck_abs_1 = np.loadtxt(directory+'/opac_planck_abs'+label+'.dat')
planck_sca_1 = np.loadtxt(directory+'/opac_planck_sca'+label+'.dat')
label = '_0.02'
planck_abs_2 = np.loadtxt(directory+'/opac_planck_abs'+label+'.dat')
planck_sca_2 = np.loadtxt(directory+'/opac_planck_sca'+label+'.dat')

plt.contour(temp, sizes, planck_abs_1, alpha=0.5)
plt.contour(temp, sizes, planck_abs_2, linestyles=':')
plt.show()

plt.contour(temp, sizes, planck_sca_1, alpha=0.5)
plt.contour(temp, sizes, planck_sca_2, linestyles=':')
plt.show()