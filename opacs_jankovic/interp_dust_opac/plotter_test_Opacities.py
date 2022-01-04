import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}
matplotlib.rc('font', **font)

###################################################
#   Plot test 1
########################

plt.clf()
plt.cla()

data = []
data.append(np.loadtxt("fayalite_F01/opac_planck_abs.dat"))
data.append(np.loadtxt("fayalite_F01/opac_planck_sca.dat"))
data.append(np.loadtxt("test_1a.dat"))
data.append(np.loadtxt("test_1b.dat"))

for d in data:
    plt.clf()
    plt.cla()
    cp = plt.contourf(d)
    plt.show()