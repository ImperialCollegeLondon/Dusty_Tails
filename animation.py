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
    T, A, X, Y, Z = line.split(',')
    time.append(float(T))
    semi.append(float(A))    
    x_pos.append(float(X))
    y_pos.append(float(Y))
    z_pos.append(float(Z))
    
    
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

#ffmpeg -start_number 0 -i fig%d.png -vcodec mpeg4 vid_out.mp4 (to get video from images, run this on command line)
