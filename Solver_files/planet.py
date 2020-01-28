import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from matplotlib.patches import Circle

matplotlib.use('Agg')

plt.ioff()

df = pd.read_csv('output.txt', sep = ",", names = ["time", "particle_id", \
                                                  "x", "y", "z"])
p_yprime = []
p_z = []


t_0 = 2.0

theta = []

for t in df['time']:
    angle = 2.0*math.pi * (t - t_0)

    theta.append(angle)

df['angles'] = theta

for angle in df['angles']:
    y_p = - math.sin(angle)*1.0
    p_yprime.append(y_p)
    p_z.append(0.0)

df['y_planet'] = p_yprime
df['z_planet'] = p_z
y_prime = []

for i in range(df['x'].size):
    yp = - math.sin(df['angles'][i])*df['x'][i] + \
    math.cos(df['angles'][i])*df['y'][i]
    y_prime.append(yp)

df['yprime'] = y_prime

i = 0

for t in df['time'].unique():
     plot_df = df[df.time == t]
     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_xlim(-1.2, 1.2)
     ax.set_ylim(-1.2, 1.2)
     ax.set_aspect('equal')
     """
     if plot_df['angles'][plot_df.index[0]] <  and \
     plot_df['angles'][plot_df.index[0]] > :


            #dust1 = plt.scatter(plot_df['yprime'], plot_df['z'], s= 4.0, c='r', alpha=0.5)
            #planet = plt.plot(plot_df['y_planet'], plot_df['z_planet'], marker = 'o', c= 'b')
            #star = plt.scatter(0.0, 0.0, s=500.0, c='orange')
            p = plt.Circle((0, 0), 0.2, color = 'orange')
            ax.add_artist(p)

            plt.savefig("fig{0:01}.png".format(i))

            plt.close()
     else:
     """
     

     dust1 = plt.scatter(plot_df['yprime'], plot_df['z'], s= 4.0, c='r', alpha=0.5)
     planet = plt.plot(plot_df['y_planet'], plot_df['z_planet'], marker = 'o', c= 'b')
     star = plt.scatter(0.0, 0.0, s=500.0, c='orange')


     plt.savefig("fig{0:01}.png".format(i))

     plt.close()

     i +=1

plt.close('all')


"""
t_0 = 2.0

theta = []

for t in time:
    angle = 2.0*math.pi * (t - t_0)

    theta.append(angle)

y_prime = []

for i in range(len(x)):
    yp = - math.sin(theta[i])*x[i] + math.cos(theta[i])*y[i]
    y_prime.append(yp)


for i in range (0, 50):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')

    planet = plt.plot([y_prime[i]], [z[i]], marker = 'o', color = 'r')

    plt.savefig("fig{0:01}.png".format(i))

    plt.close()


plt.close('all')
"""
