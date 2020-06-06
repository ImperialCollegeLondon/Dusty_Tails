import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import scipy.stats as stats
import math
from scipy.special import erf, erfc
from scipy import optimize
from textwrap import wrap
from mpl_toolkits.mplot3d import Axes3D

# cmap = plt.cm.get_cmap('magma_r')

plt.ioff()

def time(n):
    while n < 2.1:
        n = np.round(n,1)
        yield n
        n = np.round(n,1)
        n += 0.1
        n = np.round(n,1)

NR=200
NP=200
NT=200

for X in time(2.0):

    opticaldepth = open('Text_files_200_035/optical_depth_TIME=%s.txt' % X,'r')
    # R = open('Text_files_200_040/R.txt', 'r')
    phi = open('Text_files_200_035/phi.txt', 'r')
    theta = open('Text_files_200_035/theta.txt', 'r')
    # density = open('Text_files_200_040/density_TIME=%s.txt' % X, 'r')
    # rPositions = open('Text_files_200_040/rPositionsAtTime_TIME=%s.txt' % X, 'r')
    pPositions = open('Text_files_200_035/pPositionsAtTime_TIME=%s.txt' % X, 'r')
    tPositions = open('Text_files_200_035/tPositionsAtTime_TIME=%s.txt' % X, 'r')

# The axes
    # flat_R = np.loadtxt(R)
    flat_phi = np.loadtxt(phi)
    flat_theta = np.loadtxt(theta)
    # costheta = np.cos(flat_theta)

# The particle coordinates
    # flat_r = np.loadtxt(rPositions)
    flat_p = np.loadtxt(pPositions)
    flat_t = np.loadtxt(tPositions)
    cost = np.cos(flat_t)

# # Cartesian axes
    # x = flat_R*np.sin(flat_theta)*np.sin(flat_phi)
    # y = flat_R*np.sin(flat_theta)*np.cos(flat_phi)
    # z = flat_R*np.cos(flat_phi)

# # Cartesian particle coordinates
    # x_part = flat_r*np.sin(flat_t)*np.sin(flat_p)
    # y_part = flat_r*np.sin(flat_t)*np.cos(flat_p)
    # z_part = flat_r*np.cos(flat_p)

# Optical depth and density
    t = np.loadtxt(opticaldepth)
    # d = np.loadtxt(density)
    # print(X, np.amax(t))

    t3d = np.reshape(t, (NR+1,NP+1,NT+1),order='C')
    # d3d = np.reshape(d, (NR+1,NP+1,NT+1),order='C')
    # print(X, np.amax(t))


# # OPTICAL DEPTH WITH PARTICLES - spherical coordinates
    # plt.figure()
    # plt.contourf(flat_phi, flat_theta, t3d[:,:,-1])
    # # plt.colorbar()
    # plt.scatter(flat_p,flat_t,s=1,c='yellow',marker='.')
    # plt.ylim(1.55,1.6)
    # plt.xlim(-0.3,0.009)
    # plt.xlabel("Phi")
    # plt.ylabel("Theta")
    # plt.gca().set_aspect('equal')
    # plt.savefig("Optical depth (with particles) %s.png" %X)
    # plt.show()
    #
    # fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    # thing = ax1.contourf(flat_phi, flat_theta, t3d[:,:,-1])
    # fig.colorbar(thing, ax=[ax1,ax2], shrink = 0.8)
    # ax1.set_ylim(1.55,1.6)
    # ax1.set_ylabel("Theta")
    # ax1.set_aspect('equal')
    # ax2.contourf(flat_phi, flat_theta, t3d[:,:,-1])
    # ax2.scatter(flat_p,flat_t,s=1,c='yellow',marker='.')
    # ax2.set_ylim(1.55,1.6)
    # ax2.set_ylabel("Theta")
    # ax2.set_xlim(-0.3,0.009)
    # ax2.set_xlabel("Phi")
    # ax2.set_aspect('equal')
    # plt.show()

    plt.subplot(211)
    im = plt.contourf(flat_phi, flat_theta, t3d[:,:,-1])
    plt.ylim(1.55,1.6)
    plt.ylabel("Theta")
    plt.xlim(-0.3,0.009)
    plt.xlabel("Phi")
    plt.gca().set_aspect('equal')

    plt.subplot(212, facecolor='#47176a')
    # plt.contourf(flat_phi, flat_theta, t3d[:,:,-1])
    plt.scatter(flat_p,flat_t,s=1,c='white',marker='.')
    plt.ylim(1.55,1.6)
    plt.ylabel("Theta")
    plt.xlim(-0.3,0.009)
    plt.xlabel("Phi")
    plt.gca().set_aspect('equal')
    plt.tight_layout()

    plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = plt.axes([0.85, 0.1, 0.025, 0.8])
    plt.colorbar(im,cax=cax)
    plt.show()

# OPTICAL DEPTH WITHOUT PARTICLES - spherical coordinates
    # fig, ax = plt.subplots()
    # plt.contourf(flat_phi, flat_theta, t3d[:,:,-1])
    # plt.colorbar()
    # # planet = plt.Circle((0, math.pi/2), 0.002, color='r')
    # # ax.add_artist(planet)
    # plt.ylim(1.55,1.6)
    # plt.xlim(-0.3,0.009)
    # plt.xlabel("Phi")
    # plt.ylabel("Theta")
    # plt.gca().set_aspect('equal')
    # plt.savefig("Optical depth (035) Time = %s.png" %X)
    # plt.show()

# OPTICAL DEPTH WITHOUT PARTICLES - spatial?
#     plt.figure()
#     plt.contourf(flat_phi, costheta, t3d[:,:,-1], cmap = cmap)
#     plt.colorbar()
#     plt.scatter(flat_p,cost,s=1,c='r',marker='o')
#     plt.title("Optical Depth Time %s" %X)
#     plt.ylim(1.55,1.6)
#     plt.xlim(-0.3,0.009)
#     plt.xlabel("phi")
#     plt.ylabel("theta")
#     plt.gca().set_aspect('equal')
#     plt.savefig("Optical depth %s.png" %X)
#     plt.show()

# DENSITY
    # plt.figure(x)
    # plt.contourf(flat_phi, flat_theta, d3d[:,:,28])
    # plt.colorbar()
    # plt.title("Density")
    # plt.ylim(1.55,1.62)
    # plt.xlim(-0.3,0.009)
    # plt.xlabel("phi")
    # plt.ylabel("theta")
    # plt.axes().set_aspect('equal')

# # FOR LEAH
#     plt.figure()
#     plt.contourf(flat_phi, flat_theta, t3d[:,:,-1], cmap = cmap)
#     plt.colorbar()
#     plt.title("Optical Depth Time %s" %X)
#     plt.ylim(1.55,1.6)
#     plt.xlim(-0.3,0.009)
#     plt.xlabel("phi")
#     plt.ylabel("theta")
#     plt.gca().set_aspect('equal')

    # plt.savefig("For Leah %s.png" %x)
