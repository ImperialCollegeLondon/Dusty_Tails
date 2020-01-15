import numpy as np
import matplotlib.pyplot as plt


plt.ioff()


#construct dataset to plot
time = []
semi6 = []
semi7 = []
semi8 = []
semi9 = []
semi10 = []
semi11 = []
semi12 = []


ds6 = open('planet_data6.txt','r')

for line in ds6:
    line = line.strip()
    T, A, X, Y, Z= line.split(',')
    time.append(float(T)/0.02)
    semi6.append(float(A))

ds6.close()

ds7 = open('planet_data7.txt', 'r')

for line in ds7:
    line = line.strip()
    T, A, X, Y, Z= line.split(',')
    #time.append(float(T))
    semi7.append(float(A))

ds7.close()

ds8 = open('planet_data8.txt', 'r')

for line in ds8:
    line = line.strip()
    T, A, X, Y, Z= line.split(',')
    #time.append(float(T))
    semi8.append(float(A))

ds8.close()

ds9 = open('planet_data9.txt', 'r')

for line in ds9:
    line = line.strip()
    T, A, X, Y, Z= line.split(',')
    #time.append(float(T))
    semi9.append(float(A))

ds9.close()

ds10 = open('planet_data10.txt', 'r')

for line in ds10:
    line = line.strip()
    T, A, X, Y, Z= line.split(',')
    #time.append(float(T))
    semi10.append(float(A))

ds10.close()
"""
ds11 = open('data11.txt', 'r')

for line in ds11:
    line = line.strip()
    T, A= line.split(',')
    #time.append(float(T))
    semi11.append(float(A))

ds11.close()

ds12 = open('data12.txt', 'r')

for line in ds12:
    line = line.strip()
    T, A= line.split(',')
    #time.append(float(T))
    semi12.append(float(A))

ds12.close()

"""

true_a = []
for i in range(len(time)):
    true_a.append(1.0)


fig, ax = plt.subplots()


#ax.plot(time, semi6, marker = '.', c = '#0061E8', linewidth = 0.0, label = '1e-6')
ax.plot(time, semi7, marker = '.', c = '#D20000', linewidth = 0.0, label = '1e-7')
ax.plot(time, semi8, marker = '.', c = '#C838C0', linewidth = 0.0, label = '1e-8')
ax.plot(time, semi9, marker = '.', c = '#00CFFF', linewidth = 0.0, label = '1e-9')
ax.plot(time, semi10, marker = '.', c = '#00B900', linewidth = 0.0, label = '1e-10')
#ax.plot(time, semi12, marker = '.', c = '#FF9300', linewidth = 0.0, label = '1e-12')

ax.legend()
plt.xlabel('Number of orbits of dust around planet')
plt.ylabel('Semi-major axis relative error')

plt.savefig('semimajor_planet.png')
