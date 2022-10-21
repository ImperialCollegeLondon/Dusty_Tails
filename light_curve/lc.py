import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from numba import jit
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea
from matplotlib.patches import Circle

plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['xtick.minor.size'] = 5
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

df_obs = pd.read_csv('KIC1255b_observed.csv', header=None)
orb = ['0th', '1st', '2nd', 
        '3rd'
        #, '4th' 
        #,'5th'
        ]
# orb_l = ['$\\rm{0^{th}}$', '$\\rm{1^{st}}$', '$\\rm{2^{nd}}$', '$\\rm{3^{rd}}$', '$\\rm{4^{th}}$','$\\rm{5^{th}}$']
# orb_t = ['0.0', '1.0', '2.0', '3.0', '4.0', '5.0']
# t0s = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
# directory = '../simulations/Olivine_SL/KIC1255b/'
orb_l = ['$\\rm{0^{th}}$', '$\\rm{1^{st}}$', '$\\rm{2^{nd}}$', '$\\rm{3^{rd}}$', '$\\rm{4^{th}}$','$\\rm{5^{th}}$']
orb_t = ['0.0', '1.0', '2.0', 
        '3.0'
        #, '4.0'
        #, '5.0'
        ]
t0s = [0.5, 1.5, 2.5
        , 3.5
        #, 4.5
       #, 5.5
        ]
directory = '../simulations/Al2O3/KIC1255b/'

lc_file = 'light_curve.bin'
df = pd.read_csv('../input_grid.csv', header=None)
files = []
labels = []
cm = 1/2.54


plt.xlabel(r'Orbital phase')
plt.ylabel(r'Transit depth (%)')
cs = ['tab:blue', 'tab:green', 'tab:red', 'tab:purple', 'tab:pink', 'tab:orange']
i = 0
for s in df[0]:
    mdot = df[1][i]
    geom = df[2][i]
    if geom == 0:
        g = 'day'
    else:
        g = 'sph'
    n=0
    plt.figure(figsize=(26.0*cm,20.0*cm))
    plt.xlabel('Orbital phase')
    plt.ylabel('Transit depth (%)')
    plt.xlim(-0.15,0.30)
    plt.minorticks_on()
    for o in orb:
        file = directory+o+'_orb/s'+str(s)+'_mdot'+str(mdot)+'_'+g+'_t'+orb_t[n]+'/'+lc_file
        print(file)
        data = np.fromfile(file, dt)
        df_plot = pd.DataFrame(data)
        if (o=='0th'):
            y_line =[]
            for t in df_plot['time']:
                y_line.append(0.0)
            plt.plot(df_plot['time']-t0s[n], y_line,'--', linewidth=0.6, c='tab:gray', alpha=1.0)
            plt.plot(df_obs[0], -1.0*(df_obs[1]-1.0)*100, 'x', label='$\t{Kepler}$ data', ms=8.0, c='black') 
   
        plt.plot(df_plot['time']-t0s[n], -1.0*(df_plot['transit_depth']-1.0)*100,'-',label = orb_l[n]+' orbit', linewidth=1.5, ms=8.0, alpha=1.0, c = cs[n])
        n+=1
    plt.legend(loc='lower right')
    plt.gca().invert_yaxis()
    plt.savefig('../plots/Al2O3/orbits_'+str(mdot)+'mdot_'+str(s)+'micron_'+g+'.png')
    plt.close()
    i+=1


# files = [directory+'125/'+lc_file,
#         directory+'250/'+lc_file,
#         directory+'500/'+lc_file]

# files = [directory+'250_1.25e-3_1orb/'+lc_file,
#          directory+'250_2.5e-3_1orb/'+lc_file,
#          directory+'5.0e-3_250_1orb/'+lc_file]

# base = directory+'250_2.5e-3_1orb/'+lc_file
# base_data = np.fromfile(base, dt)
# df_base = pd.DataFrame(base_data)
# base = directory+'250/'+lc_file
# base_data = np.fromfile(base, dt)
# df_base = pd.DataFrame(base_data)
# cm = 1/2.54
# plt.figure(figsize=(26.0*cm,20.0*cm))
# plt.xlabel('Orbital phase')
# plt.ylabel('Transit difference')
# plt.xlim(-0.10,0.18)
# #plt.ylim(0.99, 1.001)
# plt.minorticks_on()

# plt.xlabel(r'Orbital phase')
# plt.ylabel(r'Transit depth (%)')

# #plt.min(np.arange(0.9950, 1.0006, 0.0002))

# #t0 = 3.5
# #labels = ['$\\rm{\delta_{125,250}}$','$\\rm{\delta_{500,250}}$','500', '750']
# #labels = ['1.25e-3','2.50e-3','5.00e-3']
# #markers = ['*', '^', 'x--', 'o', '--+']
# cs = ['tab:blue', 'tab:green', 'tab:red', 'tab:purple', 'tab:pink']
# i = 0
# for f in files:
#     data = np.fromfile(f, dt)
#     df = pd.DataFrame(data)
#     if (i==0):
#         y_line =[]
#         for t in df['time']:
#             y_line.append(0.0)
#         plt.plot(df['time']-t0s[i], y_line,'--', linewidth=0.6, c='black', alpha=1.0)
   
#     plt.plot(df['time']-t0s[i], -1.0*(df['transit_depth']-1.0)*100,'-',label = orb_l[i]+' orbit', linewidth=1.5, ms=8.0, alpha=1.0, c = cs[i])
#     #plt.plot(df['time']-t0, df['transit_depth'],'-',label = 'total', linewidth=2.5, ms=8.0, alpha=1.0, c = 'tab:blue')
#     #plt.plot(df['time']-t0, df['extinction'],'--',label = 'extinction', linewidth=1.5, ms=8.0, alpha=1.0, c = 'tab:red')
#     #plt.plot(df['time']-t0, df['scattering']+1.0,'-.',label = 'scattering', linewidth=1.5, ms=8.0, alpha=1.0, c = 'tab:green')
#     plt.legend(loc='lower right')
#     i+=1
# #plt.title('$\\rm{Mg_{0.8}Fe_{1.2}SiO_{4}}$, $\\rm{10~M_{\oplus}/Gyr}$, $\\rm{2~\mu m}$, spherical outflow')
# plt.gca().invert_yaxis()
# #plt.yscale('log')
# plt.savefig("../plots/Mg08Fe12SiO4/orbits_1mdot_2.5micron_sph.png")

