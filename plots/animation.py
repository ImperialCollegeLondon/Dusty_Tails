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

x_pos1 = []
y_pos1 = []
z_pos1 = []

dataset = open('data1.txt','r')

for line in dataset:
    line = line.strip()
    T, X, Y, Z = line.split(',')
    time.append(float(T))
    #semi.append(float(S))
    x_pos.append(float(X))
    y_pos.append(float(Y))
    z_pos.append(float(Z))


dataset.close()

dataset1 = open('data2.txt','r')

for line in dataset1:
    line = line.strip()
    T, X, Y, Z = line.split(',')
    #time.append(float(T))
    #semi.append(float(S))
    x_pos1.append(float(X))
    y_pos1.append(float(Y))
    z_pos1.append(float(Z))


dataset1.close()


for i in range(0, 300):
       fig = plt.figure()
       ax = fig.add_subplot(111)
       ax.set_xlim(-1.5, 1.5)
       ax.set_ylim(-1.5, 1.5)
       ax.set_aspect('equal')

       star = plt.plot([0], [0] ,marker = 'o', color='y')
       dust1 = plt.plot([x_pos[i]], [y_pos[i]],  marker='.', color='r')
       planet = plt.plot([1], [0] ,marker = 'o', color='b')
       dust2 = plt.plot([x_pos1[i]], [y_pos1[i]],  marker='.', color='g')

       plt.savefig("fig{0:01}.png".format(i))

       plt.close()


plt.close('all')

#ffmpeg -start_number 0 -i fig%d.png -vcodec mpeg4 vid_out.mp4 (to get video from images, run this on command line)
