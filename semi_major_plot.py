import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 


plt.ioff()


#construct dataset to plot
time = []
a = []

dataset = open('rk_data.txt','r')

for line in dataset:
    line = line.strip()
    X, Y = line.split(',')
    time.append(float(X))
    a.append(float(Y))
    
dataset.close()

true_a = []
for i in range(len(time)):
    true_a.append(1.0)
    

fig, ax = plt.subplots()

ax.plot(time, a, marker = '.', c = 'blue')
ax.plot(time, true_a, 'r-')

plt.savefig('semimajor.png')        

