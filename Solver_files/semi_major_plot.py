import numpy as np
import matplotlib.pyplot as plt


plt.ioff()


#construct dataset to plot
temperature = []
time = []


ds = open('sublimation_data.txt','r')

for line in ds:
    line = line.strip()
    T, t = line.split(',')
    temperature.append(float(T))
    time.append(float(t))

ds.close()


fig, ax = plt.subplots()


#ax.plot(time, semi6, marker = '.', c = '#0061E8', linewidth = 0.0, label = '1e-6')
ax.plot(temperature, time, marker = '.', c = '#D20000', linewidth = 0.0)
#ax.plot(time, semi8, marker = '.', c = '#C838C0', linewidth = 0.0, label = '1e-8')
#ax.plot(time, semi9, marker = '.', c = '#00CFFF', linewidth = 0.0, label = '1e-9')
#ax.plot(time, semi10, marker = '.', c = '#00B900', linewidth = 0.0, label = '1e-10')
#ax.plot(time, semi12, marker = '.', c = '#FF9300', linewidth = 0.0, label = '1e-12')

ax.legend()
#plt.xlabel('Number of orbits of dust around planet')
#plt.ylabel('Semi-major axis relative error')

plt.savefig('sublimation_rate.png')
