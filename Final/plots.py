import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
from scipy.special import erf, erfc
from scipy import optimize
from textwrap import wrap

plt.ioff()

def counting(n):
    while n < 10001:
        yield n
        n += 1000

def average(x): #finds the average of a list
    return sum(x)/len(x)

error = []
nop = []

NR=20
NP=20
NT=20

for x in counting(1000):
    nop.append(x)

    opticaldepth_ana = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/opticaldepth_ana_nop=%s.txt' % x,'r')
    opticaldepth_num = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/opticaldepth_num_nop=%s.txt' % x,'r')
    phi = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/phi.txt', 'r')
    theta = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/theta.txt', 'r')
    density_ana = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/density_ana_nop=%s.txt' % x, 'r')
    density_num = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/density_num_nop=%s.txt' % x, 'r')

    flat_phi = np.loadtxt(phi)
    flat_theta = np.loadtxt(theta)

    t_ana = np.loadtxt(opticaldepth_ana)
    t_num = np.loadtxt(opticaldepth_num)
    d_ana = np.loadtxt(density_ana)
    d_num = np.loadtxt(density_num)
    # print (np.size(d_num))
    # print(np.size(d_ana))
    # print (np.size(t_num))
    # print(np.size(t_ana))

    ta_norm = t_ana / max(t_ana)
    tn_norm = t_num / max(t_num)
    da = d_ana / max(d_ana)
    dn = d_num / max(d_num)

    ta3d = np.reshape(t_ana,(NR+1,NP+1,NT+1),order='C')
    tn3d = np.reshape(t_num,(NR+1,NP+1,NT+1),order='C')
    da3d = np.reshape(d_ana,(NR+1,NP+1,NT+1),order='C')
    dn3d = np.reshape(d_num,(NR,NP,NT),order='C')

    rms = []
    ana = np.asarray(t_ana)
    num = np.asarray(t_num)
    for i in range(np.size(ana)):
        if num[i] == 0.:
            rms_elem = 0.
            rms.append(rms_elem)
        else:
            e = (num[i]-ana[i])/num[i]
            rms_elem = np.sqrt(e**2.)
            rms.append(rms_elem)
    err = average(rms)
    print (err)
    error.append(float(err))

    plt.figure(x)
    plt.subplot(222)
    plt.contourf(flat_phi, flat_theta, ta3d[:,:,-1])
    plt.colorbar()
    plt.title("Analytical optical depth")
    plt.xlabel("phi")
    plt.ylabel("theta")
    plt.subplot(224)
    plt.contourf(flat_phi, flat_theta, tn3d[:,:,-1])
    plt.colorbar()
    plt.title("Numerical optical depth")
    plt.xlabel("phi")
    plt.ylabel("theta")
    plt.subplot(221)
    plt.contourf(flat_phi, flat_theta, da3d[:,:,-1])
    plt.colorbar()
    plt.title("Analytical density")
    plt.xlabel("phi")
    plt.ylabel("theta")
    plt.subplot(223)
    plt.contourf(dn3d[:,:,-1])
    plt.colorbar()
    plt.title("Numerical density")
    plt.xlabel("phi")
    plt.ylabel("theta")
    plt.show()

plt.figure()
plt.scatter(nop,error,c='r',marker='o')
plt.savefig("error.png")
plt.show()
