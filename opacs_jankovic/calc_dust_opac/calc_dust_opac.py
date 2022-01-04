import numpy as np
import radmc3dPy
import datetime

now = datetime.datetime.now()
print(now.isoformat())

directory = 'fayalite_F01'
grain_density = 4.39
#directory = 'enstatite_J98_J94_D95'
#grain_density = 3.20
#directory = 'olivine_F01'
#grain_density = 3.27
#directory = 'quartz_Z13'
#grain_density = 2.60
#directory = 'corundum_K95'
#grain_density = 4.00
#directory = 'silicon_carbide_L93'
#grain_density = 3.22
#directory = 'graphite_D84'
#grain_density = 2.16

logawidth = 0.02
str_logawidth = '_0.02'

opt_const_file = directory+'/opt_const.dat'
grain_size_grid = np.logspace(-5.0, -3.0, num=300, endpoint=True)
wavelength_grid = np.logspace(-5.0, 1.0, num=4000, endpoint=True)
'''
opt_const_data = np.loadtxt(opt_const_file)
wavelength_grid_um = opt_const_data[:,0]
wavelength_grid = wavelength_grid_um*1e-4
min_wav = 1e-5
j = 0
while wavelength_grid[j] < min_wav:
    j += 1
wavelength_grid = wavelength_grid[j:]
if wavelength_grid[0] > min_wav:
    extension = np.arange(min_wav, wavelength_grid[0], wavelength_grid[1]-wavelength_grid[0])
    wavelength_grid = np.concatenate((extension, wavelength_grid))
'''

np.savetxt(directory+'/opac_wavelength.dat', wavelength_grid)
np.savetxt(directory+'/opac_sizes.dat', grain_size_grid)

opt_const_data = np.loadtxt(opt_const_file)
output_file_abs = open(directory+'/opac_mie_abs'+str_logawidth+'.dat', 'wb')
output_file_sca = open(directory+'/opac_mie_sca'+str_logawidth+'.dat', 'wb')
output_file_gsc = open(directory+'/opac_mie_gsc'+str_logawidth+'.dat', 'wb')
for i in range(0, len(grain_size_grid)):
    grain_size = grain_size_grid[i]
    opacities = radmc3dPy.miescat.compute_opac_mie(
        fname=opt_const_file,
        matdens=grain_density,
        agraincm=grain_size,
        lamcm=wavelength_grid,
        logawidth=logawidth,
        na=40,
		extrapolate=True)	

    output_data = opacities['kabs']
    for k in output_data:
        output_file_abs.write("%e " % k)
    output_file_abs.write("\n")

    output_data = opacities['kscat']
    for k in output_data:
        output_file_sca.write("%e " % k)
    output_file_sca.write("\n")
	
    output_data = opacities['gscat']
    for k in output_data:
        output_file_gsc.write("%e " % k)
    output_file_gsc.write("\n")

    print(i)
    now = datetime.datetime.now()
    print(now.isoformat())

output_file_abs.close()
output_file_sca.close()
output_file_gsc.close()

