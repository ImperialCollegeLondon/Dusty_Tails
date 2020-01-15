import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')

plt.ioff()
#construct dataset to plot

time = []
semi = []
x_pos = []
y_pos = []
z_pos = []

dataset = open('data.txt','r')

for line in dataset:
    line = line.strip()
    T, X, Y, Z = line.split(',')
    time.append(float(T))
    #semi.append(float(S))
    x_pos.append(float(X))
    y_pos.append(float(Y))
    z_pos.append(float(Z))


dataset.close()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.set_aspect('equal')
plt.xlabel('x position in dimensionless units')
plt.ylabel('y position in dimensionless units')

plt.plot(x_pos, y_pos, '-r')

plt.savefig("total_orbit.png")
