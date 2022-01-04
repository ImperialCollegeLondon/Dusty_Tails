import numpy as np
import matplotlib.pyplot as plt

# constants
clight = 2.99792458e10
planckh = 6.6260755e-27
kB = 1.38064852e-16
sigmaSB = 5.67051e-5

# Planck function
def B_func(freq, temp):
    if planckh*freq/(kB*temp) < 1:
        exp = np.exp(planckh*freq/(kB*temp))
        return 2*planckh*freq**3/clight**2 * 1/(exp-1)
    else:
        exp = np.exp(-planckh*freq/(kB*temp))
        return 2*planckh*freq**3/clight**2 * exp/(1-exp)
	
def calc_planck_opac(lambda_grid, opac_grid, temp):
    # frequency grid and Planck function values
    freq_grid = [clight/lam  for lam in lambda_grid]
    B_grid = [B_func(freq, temp) for freq in freq_grid]
	
    # integrate planck mean opacity
    integ_B = sigmaSB*temp**4/np.pi
    integ_Bk = 0
    for i in range(1, len(lambda_grid)):
        integ_Bk += 0.5 * (B_grid[i]*opac_grid[i]+B_grid[i-1]*opac_grid[i-1])*(freq_grid[i-1]-freq_grid[i])
    planck_opac = integ_Bk/integ_B
	
    # integrate planck mean opacity with half the step to check error
    integ_Bk = 0
    for i in range(2, len(lambda_grid), 2):
        integ_Bk += 0.5 * (B_grid[i]*opac_grid[i]+B_grid[i-2]*opac_grid[i-2])*(freq_grid[i-2]-freq_grid[i])
    error = abs(planck_opac-integ_Bk/integ_B)/planck_opac
    
    if error > 2e-2:
        print("Warning! Error in calc_planck_opac is %e" % error)
        plt.loglog(lambda_grid, opac_grid)
        plt.show()
        import sys
        sys.exit()
	
    return planck_opac

# read input
directory = 'corundum_K95'
label = '_0.02'
lambda_grid = np.loadtxt(directory+'/opac_wavelength.dat')
mie_abs_grid_all = np.loadtxt(directory+'/opac_mie_abs'+label+'.dat')
mie_sca_grid_all = np.loadtxt(directory+'/opac_mie_sca'+label+'.dat')

# set temperature grid and stellar temperatures
temp_stellar = [4550, 3830]
temp_grid = np.logspace(np.log10(100), np.log10(2000), 210)
np.savetxt(directory+'/opac_temp.dat', temp_grid)

# calculate planck mean absorption
planck_grid_all = []
for i, mie_abs_grid in zip(range(len(mie_abs_grid_all)), mie_abs_grid_all):
    planck_grid = [calc_planck_opac(lambda_grid, mie_abs_grid, temp) for temp in temp_grid]
    planck_grid_all.append(planck_grid)
    print( '%.1f' % ((i+1)/len(mie_abs_grid_all)*0.5*100), '%' )
np.savetxt(directory+'/opac_planck_abs'+label+'.dat', planck_grid_all)

for temp in temp_stellar:
    planck_grid_all = []
    for mie_abs_grid in mie_abs_grid_all:
        planck_grid_all.append(calc_planck_opac(lambda_grid, mie_abs_grid, temp))
    np.savetxt(directory+'/opac_planck_abs_stellar_%.0f' % temp +label+'.dat', planck_grid_all)

# calculate planck mean scattering
planck_grid_all = []
for i, mie_sca_grid in zip(range(len(mie_sca_grid_all)), mie_sca_grid_all):
    planck_grid = [calc_planck_opac(lambda_grid, mie_sca_grid, temp) for temp in temp_grid]
    planck_grid_all.append(planck_grid)
    print( '%.1f' % ((i+1+len(mie_abs_grid_all))/len(mie_abs_grid_all)*0.5*100), '%' )
np.savetxt(directory+'/opac_planck_sca'+label+'.dat', planck_grid_all)

for temp in temp_stellar:
    planck_grid_all = []
    for mie_sca_grid in mie_sca_grid_all:
        planck_grid_all.append(calc_planck_opac(lambda_grid, mie_sca_grid, temp))
    np.savetxt(directory+'/opac_planck_sca_stellar_%.0f' % temp +label+'.dat', planck_grid_all)
