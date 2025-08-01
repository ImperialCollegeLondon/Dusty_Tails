import matplotlib.pyplot as plt
plt.rcParams["font.family"] = 'monospace'
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
plt.rcParams["font.size"] = "25"
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
orb = ['0th', 
        '1st', 
        #'2nd', 
        #'3rd'
        #, '4th' 
        #,'5th'
        ]
grid_set = [ #[120,50,350],
             [200,30,400],
             #[300,50,500]
            ]
grids = ['R=120, $\\theta$ =50, $\\phi$=350',
        'R=200, $\\theta$ =30, $\\phi$=400',
        'R=300 cells, $\\theta$ =50 cells, $\\phi$=500 cells']
# orb_l = ['$\\rm{0^{th}}$', '$\\rm{1^{st}}$', '$\\rm{2^{nd}}$', '$\\rm{3^{rd}}$', '$\\rm{4^{th}}$','$\\rm{5^{th}}$']
# orb_t = ['0.0', '1.0', '2.0', '3.0', '4.0', '5.0']
# t0s = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
# directory = '../simulations/Olivine_SL/KIC1255b/'
orb_l = ['old', 'new','$\\rm{0^{th}}$', '$\\rm{1^{st}}$', '$\\rm{2^{nd}}$', '$\\rm{3^{rd}}$', '$\\rm{4^{th}}$','$\\rm{5^{th}}$']
orb_t = [#'0.0', 
        '1.0', 
        #'2.0', 
        #'3.0'
        # '4.0'
        #, '5.0'
        ]
t0s = [#0.5, 
        1.5, 
       # 2.5
       # , 3.5
       #  4.5
       #, 5.5
        ]
#cell_size = ['2.5e-3','2.0e-3',
#'1.0e-3'
#]
cell_size = ['$\\delta_{\\rm{2.5e-3,2.0e-3}}$','$\\delta_{\\rm{1.0e-3,2.0e-3}}$']

directory = '/lustre/astro/bmce/Dusty_Tails/simulations/Mg08Fe12SiO4/dec22/'

lc_file = 'light_curve.bin'
df = pd.read_csv('../input_grid.csv', header=None)
files = []
labels = []
cm = 1/2.54

#Number of particles plots 
nparticles = ['100','250','500']
#np_l = ['$\\delta_{100,250}$', '$\\delta_{500,250}$', '500 super-particles']
np_l = ['100 super-particles', '250 super-particles', '500 super-particles']

#Number of lc grid cells 
cells = ['60  cells', '70  cells', '80  cells', '90  cells', '100 cells']
ncells = [60,70,80,90,100]
#cells = ['$\\delta_{60,80}$', '$\\delta_{70,80}$', '$\\delta_{90,80}$', '$\\delta_{100,80}$']
#cells = [60,70,90,100]

#tau_time = ['1e-2', '5.0e-3', '2.5e-3']
tau_time = [ '5e-3','2e-3']
ps = ['125', '50']
#tau_l =  ['$\\rm{250~particles, 1/100~t_{orbit}}$', '$\\rm{125~particles, 1/200~t_{orbit}}$', '$\\rm{50~~~particles, 1/500~t_{orbit}}$']
tau_l = ['$\\delta_{125,250}$', '$\\delta_{50,250}$']
#markers = ['-.',  '--', '--X', 'o', 's']
order = [2,0,1]
plt.xlabel(r'Orbital phase')
plt.ylabel(r'Transit depth absolute difference (%)')


#cs = ['#117733', '#CC6677', '#88CCEE', 'tab:purple', 'tab:pink', 'tab:orange']
# cs = ['#e41a1c','#377eb8', '#4daf4a',
#                   '#f781bf', '#a65628', '#984ea3',
#                   '#999999', '#ff7f00', '#dede00']

cs = ['tab:green', 'tab:orange','tab:orange', 'tab:red','tab:blue', 'tab:pink']
markers = ['--*','--^','^']
i = 0
sdist = 4

if (sdist==0):
  for s in df[0]:
    mdot = df[1][i]
    geom = df[2][i]
    if geom == 0:
        g = 'day'
    else:
        g = 'sph'
    n=0
    if ((s==1.5) and (mdot==2.0) and (geom==1)):

        plt.figure(figsize=(32.0*cm,28.0*cm))
        plt.xlabel('Orbital phase')
        plt.ylabel('Transit depth absolute difference(%)')
        plt.xlim(-0.15,0.20)
        plt.ylim(1.0e-7, 1.0e-2)
        plt.minorticks_on()

        # for d in dust:
        #     directory = '/lustre/astro/bmce/Dusty_Tails/simulations/'+d+'/nov22/KIC1255b/'
        #     file = directory+orb[0]+'_orb/s'+str(s)+'_mdot'+str(mdot)+'_'+g+'_t'+orb_t[0]+'/'+lc_file
        #     print(file)
        #     data = np.fromfile(file, dt)
        #     df_plot = pd.DataFrame(data)
        #     if (n==0):
        #         y_line =[]
        #         for t in df_plot['time']:
        #             y_line.append(0.0)
        #         plt.plot(df_plot['time']-t0s[n], y_line,'--', linewidth=0.6, c='tab:gray', alpha=1.0)
        #         plt.plot(df_obs[0], -1.0*(df_obs[1]-1.0)*100, 'x', label='$\t{Kepler}$ data', ms=8.0, c='black') 
    
        #     plt.plot(df_plot['time']-t0s[0], -1.0*(df_plot['transit_depth']-1.0)*100,'-',label = dust_l[n], linewidth=1.5, ms=8.0, alpha=1.0, c = cs[n])
        #     n+=1

        for step in tau_time:
            if (n==0):
                directory = '/lustre/astro/bmce/Dusty_Tails/simulations/computational_tests/nparticles'
                file = directory+'/s'+str(s)+'_mdot'+str(mdot)+'_'+g+'_2orbits_n250_Mg08Fe12SiO4/'+lc_file
                data_base = np.fromfile(file,dt)
                df_base = pd.DataFrame(data_base)
            
            directory = '/lustre/astro/bmce/Dusty_Tails/simulations/computational_tests/throwout_time/'
            file = directory+'/s'+str(s)+'_mdot'+str(mdot)+'_'+g+'_2orbits_t'+step+'_n'+ps[n]+'_Mg08Fe12SiO4/'+lc_file
            data = np.fromfile(file, dt)
            df_plot = pd.DataFrame(data)

            depths = []
            for time in df_base['time']:
                k = 0
                for t in df_plot['time']:    
                    if (abs(t-time)<(1.0e-6)):
                     depths.append(df_plot['transit_depth'][k])
                     print(df_plot['transit_depth'][k])
                    k +=1
            # if (n==0):
            #      y_line =[]
            #      for t in df_plot['time']:
            #          y_line.append(0.0)
            #      plt.plot(df_plot['time']-t0s[n], y_line,'--', linewidth=0.6, c='tab:gray', alpha=1.0) 
    
            #plt.plot(df_plot['time']-t0s[0], -1.0*(df_plot['transit_depth']-1.0)*100,markers[n],label = tau_l[n],ms=10.0, linewidth=2.0,alpha=1.0, c = cs[n], zorder=order[n])
            plt.plot(df_base['time']-t0s[0], abs(df_base['transit_depth']-depths)*100,markers[n],label = tau_l[n], linewidth=2.0, ms=10.0, alpha=1.0, c = cs[n])
            n+=1

        plt.legend(loc='lower right')
        #plt.gca().invert_yaxis()
        plt.yscale('log')
        plt.savefig('../plots/paper_plots/computational_tests/throwout_time_error.pdf', bbox_inches='tight')
        plt.savefig('../plots/paper_plots/computational_tests/throwout_time_error.png', bbox_inches='tight')
        plt.close()
    i+=1

elif (sdist==1):

    plt.figure(figsize=(26.0*cm,20.0*cm))
    plt.xlabel('Orbital phase')
    plt.ylabel('Transit depth (%)')
    plt.xlim(-0.15,0.15)
    plt.minorticks_on()
    n = 0
    for o in orb:
        if (n==0):
            file = directory+'sdist_1.75_mdot2.0_sph_t'+orb_t[n]+'/'+lc_file
            print(file)
            data = np.fromfile(file, dt)
            df_plot = pd.DataFrame(data)
        elif (n==1):
            file = directory+'sdist_1.75_mdot2.0_sph_t'+orb_t[0]+'_TAUGRIDFIX/'+lc_file
            print(file)
            data = np.fromfile(file, dt)
            df_plot = pd.DataFrame(data)

        if (n==0):
            y_line =[]
            for t in df_plot['time']:
                y_line.append(0.0)
            plt.plot(df_plot['time']-t0s[n], y_line,'--', linewidth=0.6, c='tab:gray', alpha=1.0)
            plt.plot(df_obs[0], -1.0*(df_obs[1]-1.0)*100, 'x', label='$\t{Kepler}$ data', ms=8.0, c='black') 
   
        plt.plot(df_plot['time']-t0s[0], -1.0*(df_plot['transit_depth']-1.0)*100,'-',label = orb_l[n]+' orbit', linewidth=2.0, ms=8.0, alpha=1.0, c = cs[n])
        n+=1
    plt.legend(loc='lower right', fontsize='15')
    plt.gca().invert_yaxis()
    plt.savefig('../plots/Mg08Fe12SiO4/sdist_1.75_mdot2.0_sph_taugridfix.png')
    plt.close()

elif (sdist==4):

    plt.figure(figsize=(32.0*cm,28.0*cm))
    plt.xlabel('Orbital phase')
    #plt.ylabel('Transit depth (%)')
    plt.ylabel('Transit depth absolute difference (%)')
    plt.xlim(-0.10,0.15)
    plt.ylim(1e-5,2e-2)
    plt.minorticks_on()
    n = 0
    for csize in cell_size:
        if (n==0):
            file = directory+'sdist_mu2.0_sigma1.0_mdot2.0_sph_t1.0_d2.0e-3/'+lc_file
            data_base = np.fromfile(file, dt)
            df_base= pd.DataFrame(data_base)
        if (n==1):
            # = directory+'sdist_mu2.0_sigma1.0_mdot2.0_sph_t1.0_d2.0e-3/'+lc_file
            file = directory+'sdist_mu2.0_sigma1.0_mdot2.0_sph_t1.0_r120_t50_p350/'+lc_file
        else:
            file = directory+'sdist_mu2.0_sigma1.0_mdot2.0_sph_t1.0_d2.5e-3/'+lc_file
        
        data = np.fromfile(file, dt)
        df_plot = pd.DataFrame(data)
        # if (n==0):
        #     y_line =[]
        #     for t in df_plot['time']:
        #         y_line.append(0.0)
        #     plt.plot(df_plot['time']-t0s[n], y_line,'--', linewidth=0.6, c='tab:gray', alpha=1.0)
        #     #plt.plot(df_obs[0], -1.0*(df_obs[1]-1.0)*100, 'x', label='$\t{Kepler}$ data', ms=8.0, c='black') 


        #plt.plot(df_plot['time']-t0s[0], -1.0*(df_plot['transit_depth']-1.0)*100,markers[n],label = cell_size[n], linewidth=1.5, ms=8.0, alpha=1.0, c = cs[n])
        plt.plot(df_plot['time']-t0s[0], abs(df_plot['transit_depth']-df_base['transit_depth'])*100,markers[n],label=cell_size[n], linewidth=2.0, ms=8.0, alpha=1.0, c = cs[n])
        n+=1
    plt.legend(loc='best', fontsize='25')
    #plt.gca().invert_yaxis()
    plt.yscale('log')
    plt.savefig('../plots/Mg08Fe12SiO4/sdist_mu2.0_sigma1.0_mdot2.0_sph_GRID_RES_ERR.pdf')
    plt.savefig('../plots/Mg08Fe12SiO4/sdist_mu2.0_sigma1.0_mdot2.0_sph_GRID_RES_ERR.png')
    plt.close()


sdist = 3
i=0
if (sdist==0):
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
    plt.legend(loc='lower right', fontsize='15')
    plt.gca().invert_yaxis()
    plt.savefig('../plots/Olivine_SC/light_curves/s'+str(s)+'_mdot_'+str(mdot)+'_'+g+'.png')
    plt.close()
    i+=1
    #     if (o=='0th'):
    #         y_line =[]
    #         for t in df_plot['time']:
    #             y_line.append(0.0)
    #         plt.plot(df_plot['time']-t0s[n], y_line,'--', linewidth=0.6, c='tab:gray', alpha=1.0)
    #         plt.plot(df_obs[0], -1.0*(df_obs[1]-1.0)*100, 'x', label='$\t{Kepler}$ data', ms=8.0, c='black') 
   
    #     plt.plot(df_plot['time']-t0s[n], -1.0*(df_plot['transit_depth']-1.0)*100,'-',label = orb_l[n]+' orbit', linewidth=1.5, ms=8.0, alpha=1.0, c = cs[n])
    #     n+=1

    # plt.legend(loc='lower right')
    # plt.gca().invert_yaxis()
    # plt.savefig('../plots/Mg08Fe12SiO4/light_curves/s'+str(s)+'_mdot'+str(mdot)+'_'+g+'.png')
    # plt.close()
    # i+=1



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

