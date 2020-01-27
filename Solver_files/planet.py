import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math

matplotlib.use('Agg')

plt.ioff()

time = np.arange (0.0, 5.0, 0.1)
print(time)
dataset = open('test_planet.txt','r')

x = []
y = []
z = []

for line in dataset:
    line = line.strip()
    X, Y, Z = line.split(',')
    x.append(float(X))
    y.append(float(Y))
    z.append(float(Z))

dataset.close()

t_0 = 2.2

theta = []

for t in time:
    angle = 2.0*math.pi * (t - t_0)

    theta.append(angle)

y_prime = []

for i in range(len(x)):
    yp = - math.sin(theta[i])*x[i] + math.cos(theta[i])*y[i]
    y_prime.append(yp)


for i in range (0, 35):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')

    planet = plt.plot([y_prime[i]], [z[i]], marker = 'o', color = 'r')

    plt.savefig("fig{0:01}.png".format(i))

    plt.close()


plt.close('all')
