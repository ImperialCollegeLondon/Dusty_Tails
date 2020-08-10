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


dt = np.dtype([('radius', np.float64), ('theta', np.float64), ('phi', np.float64), ('od', np.float64)])

data = np.fromfile("ray_tracer_test.bin", dt)
df = pd.DataFrame(data)

#print(df)

radii = List()
thetas = List()
phis = List()
ods = List()


for radius in df['radius'].unique():
    radii.append(radius)

for theta in df['theta'].unique():
    thetas.append(theta)

for phi in df['phi'].unique():
    phis.append(phi)

for od in df['od']:
    ods.append(od)


xs= List()
ys = List()
zs = List()


xs, ys, zs = grid(radii, thetas, phis)


#print(len(ods))
fig = plt.figure()
ax = plt.axes(projection='3d')
#ax = fig.add_subplot(111)
ax.scatter3D(df['theta'], df['phi'], ods)
#plt.scatter(df['radius'], ods)
plt.savefig('ods_3D.png')
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
ax.scatter(x, y, c=z, s=100, edgecolor='')


plt.savefig('gauss.png')
plt.close()
"""
