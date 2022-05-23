import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from numba import jit
from mpl_toolkits import mplot3d

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


dt = np.dtype([('theta', np.float64), ('phi', np.float64), ('od', np.float64)])

data = np.fromfile("./simulations/KIC1255b_Al2O3_03micro_1mdot_day_1orb_tc_optical_depth.bin", dt)
df = pd.DataFrame(data)

#print(df)
# df = df.iloc[3125:,:]
#print(df)
#print(data)
radii = []
thetas = []
phis = []
ods = []
ods_a = []
errors = []

#print(df)
for theta in df['theta'].unique():
    thetas.append(theta)

#print(thetas)
for phi in df['phi'].unique():
    phis.append(phi)

#print(len(phis))

for od in df['od']:
    if od == 0.0:
        od = np.nan
        ods.append(math.log10(od))
    elif od == np.nan:
        print("od in nan")
        ods.append(math.log10(od))
    else:
        if (od > 1.0):
            print(od)
        ods.append(math.log10(od))

# """


# for od in df['od']:
#     if math.isnan(od):
#         ods.append(0.0)
#     else:
#         ods.append(od)
# """

# #print(len(ods))

xs= []
ys = []
zs = []


# #xs, ys, zs = grid(radii, thetas, phis)


# #print(len(ods))
fig = plt.figure()
#ax = plt.axes(projection='3d')
ax = fig.add_subplot(111)
cm = 1/2.54
plt.figure(figsize=(20.0*cm,15.0*cm))
o_reshape = np.reshape(ods, (50,190), order = 'C')

# #print(o_reshape)

# #print(type(o_reshape))
levels = np.arange(-3.0, 0.3, 0.2)
plt.contourf(phis, thetas, o_reshape, levels = levels ,cmap='jet', extend='both')


# #ax.scatter3D(df['theta'], df['phi'], df['od'])
# #plt.scatter(df['radius'], ods, s= 0.01, c = "black")
# #plt.scatter(df['radius'], ods_a, s=0.01, c = "red")

cbar = plt.colorbar()
cbar.ax.set_ylabel('log(tau)')
ax.set_xlabel("phi")
ax.set_ylabel("theta")
plt.title('Optical depth at R=1.1a.\n50 cells in $\\theta$, 190 in $\\phi$ and 120 in R.\n$\\dot{M}$=1 $M_{\\oplus}$/Gyr. $s_0$=0.3$\\mu$m.')
plt.savefig('./plots/optical_depth_test.png')

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
