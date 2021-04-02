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


dt = np.dtype([('theta', np.float64), ('phi', np.float64),('ext', np.float64), ('od', np.float64)])

data = np.fromfile("./data/ray_kic1255b_3o_035_1k_25t.bin", dt)
df = pd.DataFrame(data)

#print(df)

radii = []
thetas = []
phis = []
ods = []
ods_a = []
errors = []


for theta in df['theta'].unique():
    thetas.append(theta)

print(thetas)
for phi in df['phi'].unique():
    phis.append(phi)

for od in df['od']:
    print(od)
    if od == 0.0:
        od = np.nan
        ods.append(math.log10(od))
    else:
        ods.append(math.log10(od))


"""

for od in df['od']:
    ods.append(od)
print(len(ods))
"""
xs= []
ys = []
zs = []


#xs, ys, zs = grid(radii, thetas, phis)


#print(len(ods))
fig = plt.figure()
#ax = plt.axes(projection='3d')
ax = fig.add_subplot(111)

o_reshape = np.reshape(ods, (25,150), order = 'C')

#print(o_reshape)

#print(type(o_reshape))

plt.contourf(phis, thetas, o_reshape, 20, cmap='jet', extend='both')


#ax.scatter3D(df['theta'], df['phi'], df['od'])
#plt.scatter(df['radius'], ods, s= 0.01, c = "black")
#plt.scatter(df['radius'], ods_a, s=0.01, c = "red")

cbar = plt.colorbar()
cbar.ax.set_ylabel('log(tau)')
ax.set_xlabel("phi")
ax.set_ylabel("theta")

plt.savefig('./plots/od_25_035.png')

plt.close()

"""
fig = plt.figure()
ax = fig.add_subplot(111)

plt.scatter(df ['phi'], df['extinction'])

plt.savefig('gauss.png')
plt.close()


x = df["radius"]
y = df["phi"]
z = df["extinction"]

idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=300, edgecolor='')


plt.savefig('gauss.png')
plt.close()
"""
