import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from numba import jit
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea
from matplotlib.patches import Circle

plt.rcParams["font.size"] = "18"
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

directory = '../simulations/Mg2SiO4/KIC1255b/'
#directory = '../simulations/Mg08Fe12SiO4/noparticles/'
lc_file = 'light_curve.bin'
df = pd.read_csv('../input_grid.csv', header=None)
files = []
labels = []
i = 0
for s in df[0]:
    if (df[1][i]==5.0):
      #print(s, df[1][i], df[2][i])
      if df[2][i] == 1:
        geom = 'sph'
        files.append(directory+'s'+str(s)+'_mdot'+str(df[1][i])+'_'+geom+'_t0.0/'+lc_file)
        labels.append(str(s)+'$\\rm{\mu m}$')
        #labels.append('spherical')
    #   if df[2][i] == 0:
    #     geom = 'day'
    #     files.append(directory+'s'+str(s)+'_mdot'+str(df[1][i])+'_'+geom+'_t1.0/'+lc_file)
    #     labels.append(str(df[1][i])+'$\\rm{M_{\oplus}/Gyr}$')
    #     labels.append('day-side')
    i +=1


# files = [directory+'125/'+lc_file,
#          directory+'500/'+lc_file]

# base = directory+'250/'+lc_file
# base_data = np.fromfile(base, dt)
# df_base = pd.DataFrame(base_data)
cm = 1/2.54
plt.figure(figsize=(30.0*cm,25.0*cm))
plt.xlabel('Orbital phase')
plt.ylabel('Transit difference')
plt.xlim(-0.1,0.12)
#plt.ylim(0.99, 1.001)

plt.xlabel(r'Orbital phase')
plt.ylabel(r'Transit depth (%)')

t0 = 0.5
#labels = ['$\\rm{\delta_{125,250}}$','$\\rm{\delta_{500,250}}$','500', '750']
# labels = ['1.0e-2','2.0e-2','0.5e-2','0.67e-3','0.57e-3']
#markers = ['*--', 'x--', '^', 'o', '--+']
i = 0
for f in files:
    data = np.fromfile(f, dt)
    df = pd.DataFrame(data)
    if (i==0):
        y_line =[]
        for t in df['time']:
            y_line.append(0.0)
        plt.plot(df['time']-t0, y_line,'-', linewidth=1.0, c='black', alpha=1.0)
   
    plt.plot(df['time']-t0, -1.0*(df['transit_depth']-1.0)*100,'-',label = labels[i], linewidth=2.0, alpha=0.8)
    plt.legend()
    i+=1
plt.title('Mg2SiO4, $\\rm{5M_{\oplus}/Gyr}$, spherical outflow')
plt.gca().invert_yaxis()
#plt.yscale('log')
plt.savefig("Mg2SiO4_5mdot_sph_dustsizes.png")

