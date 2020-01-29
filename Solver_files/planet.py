import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from matplotlib.patches import Circle
import seaborn

matplotlib.use('Agg')

plt.ioff()

df = pd.read_csv('output4.txt', sep = ",", names = ["time", "particle_id", \
                                                  "x", "y", "z"])
p_yprime = []
p_xprime = []
p_z = []

t_0 = 2.0

theta = []

for t in df['time']:
    angle = 2.0*math.pi * (t - t_0)
    theta.append(angle)

df['angles'] = theta

for angle in df['angles']:
    x_p = math.cos(angle)
    y_p = math.sin(angle)
    p_yprime.append(y_p)
    p_xprime.append(x_p)
    p_z.append(0.0)

df['x_planet'] = p_xprime
df['y_planet'] = p_yprime
df['z_planet'] = p_z

x_prime = []
y_prime = []

for i in range(df['x'].size):
    xp = math.cos(df['angles'][i])*df['x'][i] - \
         math.sin(df['angles'][i])*df['y'][i]

    yp =  math.sin(df['angles'][i])*df['x'][i] + \
    math.cos(df['angles'][i])*df['y'][i]

    x_prime.append(xp)
    y_prime.append(yp)

df['xprime'] = x_prime
df['yprime'] = y_prime

i = 0

for t in df['time'].unique():
     plot_df = df[df.time == t]
     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_xlim(-0.9, 0.9)
     ax.set_ylim(-0.9, 0.9)
     ax.set_aspect('equal')
     ax.set_facecolor('black')
     ax.get_xaxis().set_visible(False)
     ax.get_yaxis().set_visible(False)

     t_hours = 0.8 * 24.0 * t
     plt.title("time: %.2f hours" % t_hours)

     y_behind = []
     z_behind = []
     y_front = []
     z_front = []

     for index in plot_df.index:
         if df['xprime'][index] < 0.0:
             y_behind.append(df['yprime'][index])
             z_behind.append(df['z'][index])
         else:
             y_front.append(df['yprime'][index])
             z_front.append(df['z'][index])

     if plot_df['x_planet'][plot_df.index[0]] <  0.0:

       dust1 = plt.scatter(y_behind, z_behind, s= 10.0, c='#008080', alpha=0.3)
       planet = plt.scatter(plot_df['y_planet'], plot_df['z_planet'], s= 400.0 , c= '#33ccff')
       star = plt.scatter(0.0, 0.0, s=3000.0, c='#ffcc00')
       dust2 = plt.scatter(y_front, z_front, s=10.0, c='#008080', alpha = 0.3)


       plt.savefig("fig{0:01}.png".format(i))

       plt.close()
     else:
       dust1 = plt.scatter(y_behind, z_behind, s= 10.0, c='#008080', alpha=0.3)
       star = plt.scatter(0.0, 0.0, s=3000.0, c='#ffcc00')
       planet = plt.scatter(plot_df['y_planet'], plot_df['z_planet'], s= 400.0 , c= '#33ccff')
       dust2 = plt.scatter(y_front, z_front, s=10.0, c='#008080', alpha = 0.3)


       plt.savefig("fig{0:01}.png".format(i))

       plt.close()



     i +=1

plt.close('all')
