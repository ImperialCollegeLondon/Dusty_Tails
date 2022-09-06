import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from numba import jit
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea
from matplotlib.patches import Circle

plt.rcParams["font.size"] = "15"
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

directory = "../simulations/Mg08Fe12SiO4/throw_time/"
lc_file = 'light_curve.bin'
files = ['../simulations/Mg08Fe12SiO4/noparticles/250/'+lc_file,
         directory +'0.02/'+lc_file,
         directory +'5e-3/'+lc_file]

cm = 1/2.54
plt.figure(figsize=(25.0*cm,20.0*cm))
plt.xlabel('orbital phase')
plt.ylabel('Normalised stellar flux')
plt.xlim(-0.1,0.10)
plt.ylim(0.9992, 1.0001)

plt.xlabel(r'orbital phase')
plt.ylabel(r'Normalised stellar flux')

t0 = 0.5
#labels = ['125','250','500', '750', '1k']
labels = ['1.0e-2','2.0e-2','0.5e-2','0.67e-3','0.57e-3']
markers = ['*', 'x', '^', '--o', '--+']
i = 0
for f in files:
    data = np.fromfile(f, dt)
    df = pd.DataFrame(data)
    if (i==0):
        y_line =[]
        for t in df['time']:
            y_line.append(1.0)
        plt.plot(df['time']-t0, y_line,'-', linewidth=1.0, c='black', alpha=1.0)
   
    plt.plot(df['time']-t0, df['transit_depth'], markers[i], linewidth=1.0, label=labels[i], alpha=0.7)
    plt.legend()
    i+=1


plt.savefig("lc_throw_time.png")

