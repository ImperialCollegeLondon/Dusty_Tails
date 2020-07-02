import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from numba import jit
from numba.typed import List

@jit(nopython = True)
def to_radius_angle(y, z, sizes):
    particles_stuff = List()
    counter = 0
    for value in y:
        radius = np.sqrt(value**2.0 + z[counter]**2.0)
        phi = np.arctan(z[counter]/value)

        if radius <= 0.20:
            particles_stuff.append([radius, phi, sizes[counter]])
        counter += 1
    return particles_stuff

@jit(nopython = True)
def patch(radii, phis, area):
    for r in range(1, len(radii)):
        for a in range(1, len(phis)):
            area[r][a] = 0.5*(phis[a]-phis[a-1])*(radii[r]**2. - radii[r-1]**2.)
    return area

@jit(nopython = True)
def where_grid(particles_stuff, radii, phis, area):
    optical_depths = np.zeros((100, 100), dtype = np.float64)
    for particle in particles_stuff:
      for r in range(1, len(radii)):
        if (radii[r-1] < particle[0] < radii[r]):
          for a in range(1, len(phis)):
            if (phis[a-1] < particle[1] < phis[a]):
                s = (np.pi * (particle[2] / (1.30*10.**(11.)))**(2.) * (8.6*10.**(24.))) / area[r][a]

                optical_depths[r][a] = optical_depths[r][a] + s
    return optical_depths

@jit(nopython = True)
def flux(area, optical_depths, radii, phis):
    f = 0.0
    for r in range(1, len(radii)):
        for a in range(1, len(phis)):
            #print(area[r][a])
            f = f + ((area[r][a] * np.exp(-1.0*optical_depths[r][a])) / (np.pi * 0.20**(2.)))

    return f


dt = np.dtype([('time', np.float64), ('id', np.int64), ('x', np.float64), \
('y', np.float64), ('z', np.float64), ('size', np.float64), ('mass', np.float64)])

data = np.fromfile("output_k222_035.bin", dt)
df = pd.DataFrame(data)

p_yprime = []
p_xprime = []
p_z = []

t_0 = 0.0

theta = []


for t in df['time']:
    angle = 2.0*math.pi * (t - t_0)
    theta.append(angle)

df['angles'] = theta

for angle in df['angles']:
    x_p = math.cos(angle) * math.sin(1.38)
    y_p = math.sin(angle)
    z_p = -math.cos(angle) * math.cos(1.38)
    p_yprime.append(y_p)
    p_xprime.append(x_p)
    p_z.append(z_p)

df['x_planet'] = p_xprime
df['y_planet'] = p_yprime
df['z_planet'] = p_z

x_prime = []
y_prime = []

for i in range(df['x'].size):
    xp = math.cos(df['angles'][i])*df['x'][i] - \
         math.sin(df['angles'][i])*df['y'][i]

    yp =  math.sin(df['angles'][i])*df['x'][i] + \
    math.cos(df['angles'][i])*df['y'][i]

    x_prime.append(xp)
    y_prime.append(yp)

df['xprime'] = x_prime
df['yprime'] = y_prime

x_dprime = []
z_dprime = []

for i in range(df['xprime'].size):
    incl = 1.38

    xdp = df['xprime'][i]* math.sin(incl) + df['z'][i]*math.cos(incl)
    zdp = -df['xprime'][i]* math.cos(incl) + df['z'][i]*math.sin(incl)
    x_dprime.append(xdp)
    z_dprime.append(zdp)

df['x_double_prime'] = x_dprime
df['z_double_prime'] = z_dprime


radii = List()
phis = List()

radii = np.linspace(0.0, 0.20, num = 100)
phis = np.linspace(0.0, 2.*np.pi, num = 100)

area = np.zeros((100, 100), dtype = np.float64)

area_new = patch(radii, phis, area)

t_flux = []

total_flux = []

for t in df['time'].unique():
  if (t > 1.75) and (t < 2.25):

    #print("this is running ")
    plot_df = df[df.time == t]

    t_hours = 9.146 * t
    #plt.title("time: %.2f hours" % t_hours)

    y_behind = []
    z_behind = []
    y_front = List()
    z_front = List()
    p_sizes = List()

    for index in plot_df.index:
         if (plot_df['x_double_prime'][index] > 0.0):

             y_front.append(plot_df['yprime'][index])
             z_front.append(plot_df['z_double_prime'][index])
             p_sizes.append(plot_df['size'][index])

    particles = List()

    particles = to_radius_angle(y_front,z_front,p_sizes)
    depths = where_grid(particles, radii, phis, area_new)

    #print(depths)
    total = flux(area, depths, radii, phis)
    print(t, total)

    t_flux.append(t-2.)
    total_flux.append(total)


#plt.close('all')
final_data = {'time': t_flux, 'normalised flux': total_flux}
final_df = pd.DataFrame(final_data, columns = ['time', 'normalised flux'])

final_df.to_csv("output_k222_035_max.txt", sep = ',')
fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(t_flux, total_flux)
ax.set_ylim(0.990, 1.001)
ax.set_xlim(1.9, 2.2)
ax.set_xlabel('Orbital Phase')
ax.set_ylabel('Normalised Stellar Flux')
ax.autoscale()


plt.savefig('lightcurve_k222_035_max.png')
plt.close()
