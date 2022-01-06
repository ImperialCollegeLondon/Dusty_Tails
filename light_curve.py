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

        if radius <= 0.24:
            if math.isnan(sizes[counter]) == False:
                #print("particle within star")
                particles_stuff.append([radius, phi, sizes[counter]])
        counter += 1
    return particles_stuff

@jit(nopython = True)
def patch(radii, phis, area):
    total_area = 0.0
    for r in range(1, len(radii)-1):
        for a in range(1, len(phis)-1):
            area[r][a] = 0.5*(phis[a+1]-phis[a])*(radii[r+1]**2. - radii[r]**2.)
            total_area = total_area + area[r][a]
    print('total_area ', total_area)
    return area

@jit(nopython = True)
def where_grid(particles_stuff, radii, phis, area):
    optical_depths = np.zeros((200, 200), dtype = np.float64)
    centre_patch = 0.0
    centre_area = np.pi*radii[1]**2.0

    for particle in particles_stuff:
      if (particle[0]<=radii[1]):
          s = (np.pi * (particle[2] / (1.929*10.**(11.)))**(2.) * (1.914*10.**(22.))) / centre_area
          centre_patch = centre_patch + s
      else:
        for r in range(1, len(radii)-1):
            if (radii[r] <= particle[0] < radii[r+1]):
                for a in range(1, len(phis)-1):
                    if (phis[a] <= particle[1] < phis[a+1]):
                        s = (np.pi * (particle[2] / (1.929*10.**(11.)))**(2.) * (1.914*10.**(22.))) / area[r][a]
                        optical_depths[r][a] = optical_depths[r][a] + s
    return optical_depths, centre_patch

@jit(nopython = True)
def flux(area, optical_depths,centre_tau, radii, phis):
    f = (np.pi*radii[1]**2.0 * np.exp(-1.0*centre_tau) )/(np.pi * 0.24**(2.))
    for r in range(1, len(radii)):
        for a in range(1, len(phis)):
            f = f + ((area[r][a] * np.exp(-1.0*optical_depths[r][a])) / (np.pi * 0.24**(2.)))

    return f


dt = np.dtype([('time', np.float64), ('id', np.int64), ('x', np.float64), \
('y', np.float64), ('z', np.float64), ('size', np.float64), ('mass', np.float64), ('tau', np.float64)])

data = np.fromfile("./data/KIC1255b_sph_03micro_3orb_3mdot_tauconst_25t.bin", dt)
df = pd.DataFrame(data)

print("Data frame has been created.")

p_yprime = []
p_xprime = []
p_z = []

t_0 = np.pi
theta = []

for t in df['time']:
    angle = 2.0*math.pi * (t - t_0)
    #print("angle", angle)
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


print("Got planet position.")
x_prime = []
y_prime = []

x_dp = []
z_dp = []

for i in range(df['x'].size):
    xp = math.cos(df['angles'][i])*df['x'][i] - \
         math.sin(df['angles'][i])*df['y'][i]

    yp =  math.sin(df['angles'][i])*df['x'][i] + \
    math.cos(df['angles'][i])*df['y'][i]

    x_prime.append(xp)
    y_prime.append(yp)

print("Got dust particles position projections.")

df['xprime'] = x_prime
df['yprime'] = y_prime

for i in range(df['xprime'].size):
    incl = 1.38

    xdp = df['xprime'][i]* math.sin(incl) + df['z'][i]*math.cos(incl)
    zdp = -df['xprime'][i]* math.cos(incl) + df['z'][i]*math.sin(incl)
    x_dp.append(xdp)
    z_dp.append(zdp)

df['x_double_prime'] = x_dp
df['z_double_prime'] = z_dp

print("Got double projections.")

radii = List()
phis = List()

radii = np.linspace(0.0, 0.24, num = 200)
phis = np.linspace(0.0, 2.*np.pi, num = 200)

area = np.zeros((200, 200), dtype = np.float64)

area_new = patch(radii, phis, area)

t_flux = []

total_flux = []

i = 0

for t in df['time'].unique():

    print("At time ", t)
    plot_df = df[df.time == t]

    t_hours = 15.68 * t
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

    if len(y_front) == 0:
        t_flux.append(t)
        total_flux.append(1.0)
        print("Flux 1.0")
    else:
        #print("particles ", len(y_front))
        particles = to_radius_angle(y_front,z_front,p_sizes)
        depths, centre_tau = where_grid(particles, radii, phis, area_new)

        total = flux(area, depths, centre_tau,radii, phis)

        print("Flux ", total)


        t_flux.append(t)
        total_flux.append(total)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_ylim(0.95, 1.002)
    ax.set_xlim(0.0, 3.0)

    plt.plot(t_flux, total_flux)



    plt.savefig("./plots/LC_TC_25_3mdot_03micro_200grid{0:01}.png".format(i))

    plt.close()

    i +=1




plt.close('all')
"""
final_data = {'time': t_flux, 'normalised flux': total_flux}
final_df = pd.DataFrame(final_data, columns = ['time', 'normalised flux'])

final_df.to_csv("output.txt", sep = ',')
fig = plt.figure()
ax = fig.add_subplot(111)

plt.plot(t_flux, total_flux)
ax.set_ylim(0.990, 1.001)
ax.set_xlim(2.75, 3.25)
ax.autoscale()


plt.savefig('lightcurve3.png')
plt.close()
"""
