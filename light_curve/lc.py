import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from numba import jit
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea
from matplotlib.patches import Circle

plt.rcParams["font.size"] = "12"
@jit(nopython = True)
def grid(r, t, p) :
    x = []
    y = []
    z = []

    for i in range(0, len(r)):
        for j in range(0, len(t)):
            for k in range(0, len(p)):
                x.append( r[i] * math.sin(t[j]) * math.cos(p[k]))
                y.append( r[i] * math.sin(t[j]) * math.sin(p[k]))
                z.append( r[i] * math.cos(t[j]))

    return x, y, z


# dt = np.dtype([('theta', np.float64), ('phi', np.float64), ('od', np.float64)])

# data = np.fromfile("./simulations/KIC1255b_Al2O3_1micro_1mdot_sph_1orb/optical_depth.bin", dt)
# df = pd.DataFrame(data)


dt = np.dtype([('time', np.float64), ('extinction', np.float64), ('scattering', np.float64), ('transit_depth', np.float64)])
dt2 = np.dtype([('lambda', np.float64), ('absorption', np.float64)])

#directory = '../simulations/MgSiO3/KIC1255b/0th_orb/'
directory = '../simulations/Mg08Fe12SiO4/noparticles/'
lc_file = 'light_curve.bin'
#df = pd.read_csv('../input_grid.csv', header=None)
files = []
labels = []
# i = 0
# for s in df[0]:
#     #print(s, df[1][i], df[2][i])
#     if df[1][i] == 2.0:
#      if df[2][i] == 1:
#         geom = 'sph'
#         files.append(directory+'s'+str(s)+'_mdot'+str(df[1][i])+'_'+geom+'_t0.0/'+lc_file)
#         labels.append(str(s)+'$\\rm{\mu m}$')
#     i +=1


files = [directory+'125/'+lc_file,
         directory+'500/'+lc_file]

base = directory+'250/'+lc_file
base_data = np.fromfile(base, dt)
df_base = pd.DataFrame(base_data)
cm = 1/2.54
plt.figure(figsize=(7.5,6))
plt.xlabel('Orbital phase')
plt.ylabel('Transit difference')
plt.xlim(-0.1,0.15)
#plt.ylim(0.99, 1.001)

plt.xlabel(r'Orbital phase')
plt.ylabel(r'Transit difference')

t0 = 0.5
labels = ['$\\rm{\delta_{125,250}}$','$\\rm{\delta_{500,250}}$','500', '750']
# labels = ['1.0e-2','2.0e-2','0.5e-2','0.67e-3','0.57e-3']
markers = ['*--', 'x--', '^', 'o', '--+']
i = 0
for f in files:
    data = np.fromfile(f, dt)
    df = pd.DataFrame(data)
    # if (i==0):
    #     y_line =[]
    #     for t in df['time']:
    #         y_line.append(1.0)
    #     plt.plot(df['time']-t0, y_line,'-', linewidth=1.0, c='black', alpha=1.0)
   
    plt.plot(df['time']-t0, abs(df['transit_depth']-df_base['transit_depth']), markers[i],label = labels[i], linewidth=1.0, alpha=0.8)
    plt.legend()
    i+=1
#plt.title('MgSiO3. $\\rm{\dot{M} = 2.0 M_{\oplus}/Gyr }$. Spherical outflow.')
plt.yscale('log')
plt.savefig("no_particles_error.pdf")

