import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
from scipy.special import erf, erfc
from scipy import optimize
from textwrap import wrap

plt.ioff()

mean = 1.
variance = (0.2*mean)**2
sd = math.sqrt(variance)
z = np.linspace(0., 2.0, 1000.)
g = 1.44e-14*stats.norm.pdf(z,mean,sd) #analytic gaussian curve
f = (1.05/4.)*erf((z-mean)/(np.sqrt(2)*sd))+(1.05/4.) #analytic erf function

def factor_two(n):
    while n < 21:
        yield n
        n += 1

NR = []
NRu = []

for i in factor_two(11):

    x = []
    y = []
    NRu.append(i)

    dataset = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/gaussian_NRu=%s.txt' % i,'r')

    for line in dataset:
        line = line.strip()
        X,Y = line.split(',')
        x.append(float(X))
        y.append(float(Y))

    dataset.close()
    plt.figure()
    plt.scatter(x,y,c='r',marker='o') #numerical density
    plt.plot(z,g) #analytical gaussian function
    plt.xlabel("Rb/a")
    plt.ylabel("Density")
    plt.xlim(0.,2.0)
    plt.ylim(0.,1.4e-14)
    plt.title("Density as a function of Rb (dimensionless)")
    plt.savefig("gaussianu_%s.png" % i)


    a = []
    b = []

    dataset = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/optical_depth_NRu=%s.txt' % i, 'r')

    for line in dataset:
        line = line.strip()
        A,B = line.split(',')
        a.append(float(A))
        b.append(float(B))

    #plotting optical depth as a function of Rb
    dataset.close()
    plt.figure()
    plt.scatter(a,b,c='r',marker='o') #numerical optical depth
    plt.plot(z,f) #analytical erf function
    plt.xlabel("Rb/a")
    plt.ylabel("Optical Depth")
    plt.title("Optical depth as a function of Rb (dimensionless)")
    plt.savefig("optical_depthu_%s.png" % i)
