import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.integrate as integrate
import math
from scipy.special import erf, erfc
from scipy import optimize

plt.ioff()

#construct dataset to plot
x = []
y = []

dataset = open('gaussian_grid.txt','r')

for line in dataset:
    line = line.strip()
    X,Y = line.split(',')
    x.append(float(X))
    y.append(float(Y))

dataset.close()

plt.figure()
plt.scatter(x,y,c='r',marker='o')
plt.savefig("example_gaussian_grid.png")

#DENSITY construct dataset to plot
a = []
b = []

dataset = open('inverse_gaussian_grid.txt','r')

for line in dataset:
    line = line.strip()
    A,B = line.split(',')
    a.append(float(A))
    b.append(float(B))

dataset.close()

plt.figure()
plt.scatter(a,b,c='r',marker='o')
plt.savefig("example_inverse_gaussian_grid.png")

#construct dataset to plot
p = []
q = []

dataset = open('DR.txt','r')

for line in dataset:
    line = line.strip()
    P,Q = line.split(',')
    p.append(float(P))
    q.append(float(Q))

dataset.close()

plt.figure()
plt.xlabel("Rb")
plt.ylabel("DR in terms of a")
plt.scatter(p,q,c='r',marker='o')
plt.savefig("example_DR.png")




z = []
f = []

dataset = open('xRa.txt','r')

for line in dataset:
    line = line.strip()
    Z,F = line.split(',')
    z.append(float(Z))
    f.append(float(F))

dataset.close()

plt.figure()
plt.xlabel("x")
plt.ylabel("Ra")
plt.scatter(z,f,c='r',marker='o')
plt.savefig("example_xRa.png")
