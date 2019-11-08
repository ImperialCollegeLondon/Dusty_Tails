import numpy as np
import matplotlib.pyplot as plt


plt.ioff()


#construct dataset to plot
Rb = []
t = []

dataset = open('optical_depth.txt','r')

for line in dataset:
    line = line.strip()
    RB,T = line.split(',')
    Rb.append(float(RB))
    t.append(float(T))

dataset.close()

plt.plot(Rb,t)
plt.xlabel("Distance from star /Rsun (m)")
plt.ylabel("Optical depth")
plt.xlim(0., 2.1)
plt.ylim(0., 1.1e+15)
plt.savefig("optical_depth.png")
