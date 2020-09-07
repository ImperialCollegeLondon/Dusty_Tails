import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from numba import jit
from numba.typed import List

@jit(nopython = True)
def spherical(x, y, z):
    particle_spherical = List()

    for i in range(0, len(x)):
        radius = np.sqrt(x[i]**2.0 + y[i]**2.0 + z[i]**2.0)
        theta = np.arccos(z[i]/radius)
        phi = np.arctan(y[i]/x[i])

        particle_spherical.append([radius, theta, phi, x[i], y[i], z[i]])

    return particle_spherical


dt = np.dtype([('time', np.float64), ('id', np.int64), ('x', np.float64), \
('y', np.float64), ('z', np.float64), ('size', np.float64), ('mass', np.float64)])

data = np.fromfile("kic_1255b_035_spherical.bin", dt)
df = pd.DataFrame(data)


i = 0

all_r_max = []
all_r_min = []

all_theta_max = []
all_theta_min = []

all_phi_max = []
all_phi_min = []

for t in df['time'].unique():

    plot_df = df[df.time == t]

    t_hours = 15.68 * t

    particles = List()
    x_now = List()
    y_now = List()
    z_now = List()

    for index in plot_df.index:
             x_now.append(plot_df['x'][index])
             y_now.append(plot_df['y'][index])
             z_now.append(plot_df['z'][index])

    particles = List()

    particles = spherical(x_now, y_now, z_now)

    radii = List()
    thetas = List()
    phis = List()

    for particle in particles:
        #print(particle)
        radii.append(particle[0])
        thetas.append(particle[1])
        #print(particle[1])
        phis.append(particle[2])

    all_r_max.append(max(radii))
    all_r_min.append(min(radii))
    all_theta_max.append(max(thetas))
    all_theta_min.append(min(thetas))
    all_phi_max.append(max(phis))
    all_phi_min.append(min(phis))


print("max radius" , max(all_r_max))
print("min radius" , min(all_r_min))
print("max theta ", max(all_theta_max))
print("min theta ", min(all_theta_min))
print("max phi ", max(all_phi_max))
print("min phi", min(all_phi_min))







plt.close('all')
