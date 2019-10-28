import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 


plt.ioff()


#construct dataset to plot
x_pos = []
y_pos = []
z_pos = []

dataset = open('test_data.txt','r')

for line in dataset:
    line = line.strip()
    X, Y, Z = line.split(',')
    x_pos.append(float(X))
    y_pos.append(float(Y))
    
dataset.close()


for i in range(0, 100):
       fig = plt.figure()
       ax = fig.add_subplot(111)
       ax.set_xlim(-1.0, 1.0)
       ax.set_ylim(-1.0, 1.0) 
       ax.set_aspect('equal')
        
       star = plt.plot([0], [0] ,marker = 'o', color='b')
       dust = plt.plot([x_pos[i]], [y_pos[i]],  marker='o', color='r')
        
       plt.savefig("fig{0:01}.png".format(i))
       
       plt.close()
        
plt.close('all')
