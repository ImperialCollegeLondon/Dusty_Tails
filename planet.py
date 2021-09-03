import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import math
import pandas as pd
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

matplotlib.use('Agg')

jet = cm.get_cmap('jet', 30)

plt.ioff()

dt = np.dtype([('time', np.float64), ('id', np.int64), ('x', np.float64), \
('y', np.float64), ('z', np.float64), ('size', np.float64), ('mass', np.float64), ('tau', np.float64)])

data = np.fromfile("./data/KIC1255b_theta08.bin", dt)
df = pd.DataFrame(data)
#print(df['x'])
#print(df['time'])
p_yprime = []
p_xprime = []
p_z = []

t_0 = 0.0

theta = []
"""
for t in df['time'].unique():
        
        if (t>1.495) and (t<1.505):
         plot_df = df[df.time == t]
         counter=0
         for tau in plot_df['tau']:
              if(tau>1.0) and (tau<20.0):
                  counter +=1
print(counter)
         #plt.hist(df['tau'])
         #plt.savefig("histogram.png")


for i in df['tau']:
    if (i>0.1):
        #print(i)

"""

#tau_norm = [float(i)/(df['tau'].sum()) for i in df['tau']]

#print(tau_norm)

#df['tau_norm'] = tau_norm


for t in df['time']:
    angle = 2.0*math.pi * (t - t_0)
    theta.append(angle)

df['angles'] = theta

for angle in df['angles']:
    x_p = math.cos(angle) * math.sin(1.38)
    y_p = math.sin(angle)
    z_p = -math.cos(angle) * math.cos(1.38)
    p_yprime.append(y_p)
    p_xprime.append(x_p)
    p_z.append(z_p)

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

x_dprime = []
z_dprime = []

for i in range(df['xprime'].size):
    incl = 1.38

    xdp = df['xprime'][i]* math.sin(incl) + df['z'][i]*math.cos(incl)
    zdp = -df['xprime'][i]* math.cos(incl) + df['z'][i]*math.sin(incl)
    x_dprime.append(xdp)
    z_dprime.append(zdp)

df['x_double_prime'] = x_dprime
df['z_double_prime'] = z_dprime


i = 0

for t in df['time'].unique():

     plot_df = df[df.time == t]
     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.set_xlim(-1.0, 1.0)
     ax.set_ylim(-1.0, 1.0)
     ax.set_aspect('equal')
     #ax.set_facecolor('black')
     ax.get_xaxis().set_visible(False)
     ax.get_yaxis().set_visible(False)

     t_hours = 15.68* t
     #plt.title("time: %.2f hours" % t_hours)

     y_behind = []
     z_behind = []
     y_front = []
     z_front = []

     for index in plot_df.index:
         if (plot_df['x_double_prime'][index] < 0.0):
             y_behind.append(plot_df['yprime'][index])
             z_behind.append(plot_df['z_double_prime'][index])
         else:
             y_front.append(plot_df['yprime'][index])
             z_front.append(plot_df['z_double_prime'][index])




     if plot_df['x_planet'][plot_df.index[0]] <  0.0:

       dust1 = plt.scatter(y_behind, z_behind, s= 0.1,  alpha=0.01, c='#008080')
       planet = plt.scatter(plot_df['y_planet'], plot_df['z_planet'], s=20.0 , c= '#C6492B')
       star = plt.scatter(0.0, 0.0, s=10000.0, c='#ffcc00')
       dust2 = plt.scatter(y_front, z_front, s=0.1, alpha = 0.01, c='#008080')


       plt.savefig("./plots/kic_08test{0:01}.png".format(i))

       plt.close()
     else:
       dust1 = plt.scatter(y_behind, z_behind, s= 0.1, alpha=0.01, c='#008080')
       star = plt.scatter(0.0, 0.0, s=10000.0, c='#ffcc00')
       planet = plt.scatter(plot_df['y_planet'], plot_df['z_planet'], s= 20.0 , c= '#C6492B')
       dust2 = plt.scatter(y_front, z_front, s=0.1, alpha = 0.01, c='#008080')


       plt.savefig("./plots/kic_08test{0:01}.png".format(i))

       plt.close()



     i +=1

plt.close('all')
