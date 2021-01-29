import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from numba import jit
from numba.typed import List
from mpl_toolkits import mplot3d

@jit(nopython = True)
def grid(r, t, p) :
    x = List()
    y = List()
    z = List()

    for i in range(0, len(r)):
        for j in range(0, len(t)):
            for k in range(0, len(p)):
                x.append( r[i] * math.sin(t[j]) * math.cos(p[k]))
                y.append( r[i] * math.sin(t[j]) * math.sin(p[k]))
                z.append( r[i] * math.cos(t[j]))

    return x, y, z


dt = np.dtype([('theta', np.float64), ('phi', np.float64),('ext', np.float64), ('od', np.float64)])

data = np.fromfile("ray_tracer_testing_1000.bin", dt)
df = pd.DataFrame(data)

#print(df)

radii = List()
thetas = List()
phis = List()
ods = List()
ods_a = List()
errors = List()


for theta in df['theta'].unique():
    thetas.append(theta)

for phi in df['phi'].unique():
    phis.append(phi)

for od in df['od']:
    if od == 0.0:
        od = np.nan
        ods.append(math.log10(od))
    else:
        ods.append(math.log10(od))

print (df)
"""

for od in df['od']:
    ods.append(od)
print(len(ods))
"""
xs= List()
ys = List()
zs = List()


#xs, ys, zs = grid(radii, thetas, phis)


#print(len(ods))
fig = plt.figure()
#ax = plt.axes(projection='3d')
ax = fig.add_subplot(111)

o_reshape = np.reshape(ods, (120,120), order = 'C')

print(type(o_reshape))

plt.contourf(phis, thetas, o_reshape, 20, cmap='jet')

plt.colorbar()
#ax.scatter3D(df['theta'], df['phi'], df['od'])
#plt.scatter(df['radius'], ods, s= 0.01, c = "black")
#plt.scatter(df['radius'], ods_a, s=0.01, c = "red")

ax.set_xlabel("theta")
ax.set_ylabel("phi")

plt.savefig('od_test.png')

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
