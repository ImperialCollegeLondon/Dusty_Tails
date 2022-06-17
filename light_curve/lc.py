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

data_1 = np.fromfile("./data/KIC1255b_Al2O3_day_1mdot_1micro_1orb_tau_lightcurve.bin", dt)
data_2 = np.fromfile("./data/KIC1255b_Al2O3_day_1mdot_1micro_2orb_tau_lightcurve.bin", dt)

df_1 = pd.DataFrame(data_1)
df_2 = pd.DataFrame(data_2)
c = 0

for time in df_2['time']:
    df_2['time'][c] = time + 1.0;
    print(df_2['scattering'][c])
    c +=1
df = pd.concat([df_1, df_2])
cm = 1/2.54
plt.figure(figsize=(20.0*cm,18.0*cm))
plt.xlabel('phase')
plt.ylabel('transit depth')
#plt.xlim(-0.5, 0.5)
y_line = []
for time in df['time']:
    # print(time)
    y_line.append(1.0)
plt.plot(df['time'], y_line, '--', linewidth = 1.0, alpha=0.7)
plt.plot(df['time'], df['transit_depth'], '-', linewidth=2.0, label='total')
plt.plot(df['time'], df['extinction'], '--', linewidth=1.5, label='extinction', alpha=0.7)
plt.plot(df['time'], df['scattering']+1.0, '--', linewidth=1.5, label='scattering', alpha=0.7)
#print(df['depth'].idxmin())
#print(df['time'][df['transit_depth'].idxmin()])
plt.legend()
plt.savefig("light_curve_KIC1255b_day_Al2O3_1micro.png")
# #print(df)
# # df = df.iloc[3125:,:]
# #print(df)
# #print(data)
# radii = []
# thetas = []
# phis = []
# ods = []
# ods_a = []
# errors = []

# #print(df)
# for theta in df['theta'].unique():
#     thetas.append((0.05/40.)*theta + 1.55)

# #print(thetas)
# for phi in df['phi'].unique():
#     phis.append((0.19/120.)*phi - 0.18)

# #print(len(phis))

# for od in df['od']:
#     if od == 0.0:
#         od = np.nan
#         ods.append(math.log10(od))
#     elif od == np.nan:
#         print("od in nan")
#         ods.append(math.log10(od))
#     else:
#         if (od > 1.0):
#             print(od)
#         ods.append(math.log10(od))

# # """


# # for od in df['od']:
# #     if math.isnan(od):
# #         ods.append(0.0)
# #     else:
# #         ods.append(od)
# # """

# # #print(len(ods))

# xs= []
# ys = []
# zs = []


# # #xs, ys, zs = grid(radii, thetas, phis)


# # #print(len(ods))
# fig = plt.figure()
# #ax = plt.axes(projection='3d')
# ax = fig.add_subplot(111)
# cm = 1/2.54
# plt.figure(figsize=(30.0*cm,20.0*cm))
# o_reshape = np.reshape(ods, (40,120), order = 'C')

# # #print(o_reshape)

# # #print(type(o_reshape))
# levels = np.arange(-3.0, 0.3, 0.2)
# plt.contourf(phis, thetas, o_reshape, levels = levels ,cmap='jet', extend='both')


# # #ax.scatter3D(df['theta'], df['phi'], df['od'])
# # #plt.scatter(df['radius'], ods, s= 0.01, c = "black")
# # #plt.scatter(df['radius'], ods_a, s=0.01, c = "red")

# cbar = plt.colorbar()
# cbar.ax.set_ylabel('log($\\tau$)')

# p = Circle((0.0, 1.57), 20)
# ax.add_patch(p)

# plt.plot(0.0, math.pi/2., 'x', alpha=1.0, markersize=10, color='black')

# plt.title('Optical depth at R=1.1$a_{p}$. 40 cells in $\\theta$, 152 in $\\phi$ and 96 in R.\n$\\dot{M}$=1 $M_{\\oplus}$/Gyr. $s_0$=1.0$\\mu$m.\n The ranges for this are: 0.98<R<1.10, 1.55<theta<1.60, -0.18<phi<0.01 \n 500 super-particles were thrown out every 100th of an orbit.', 
# fontsize = '10')
# plt.xlabel("$\\phi$ (rad)")
# plt.ylabel("$\\theta$ (rad)")
# plt.savefig('./plots/tau_1micro_Al2O3_1mdot_sph.png')

# plt.close()

# """
# fig = plt.figure()
# ax = fig.add_subplot(111)

# plt.scatter(df ['phi'], df['extinction'])

# plt.savefig('gauss.png')
# plt.close()


# x = df["radius"]
# y = df["phi"]
# z = df["extinction"]

# idx = z.argsort()
# x, y, z = x[idx], y[idx], z[idx]

# fig, ax = plt.subplots()
# ax.scatter(x, y, c=z, s=300, edgecolor='')


# plt.savefig('gauss.png')
# plt.close()
# """
