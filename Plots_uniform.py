import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
from scipy.special import erf, erfc
from scipy import optimize

plt.ioff()

#defining variables
mean = 1.
variance = (0.2*mean)**2
sd = math.sqrt(variance)
z = np.linspace(0., 2.0, 1000.)
g = 40.3877*stats.norm.pdf(z,mean,sd) #analytic gaussian curve
f = (1.75668e+15/2.)*erf((z-mean)/(np.sqrt(2)*sd))+(1.75668e+15/2.) #analytic erf function

#defining functions
def average(x): #finds the average of a list
    return sum(x)/len(x)

def factor_two(n):
    while n < 1000:
        yield n
        n += 1

#defining lists to plot error against NR
error = []
NR = []

#for loop to produce optical depth plots for different NR
for i in factor_two(1):

    #DENSITY construct dataset to plot
    x = []
    y = []
    NR.append(i)

    dataset = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/gaussian_NR=%s.txt' % i,'r')

    for line in dataset:
        line = line.strip()
        X,Y = line.split(',')
        x.append(float(X))
        y.append(float(Y))

    dataset.close()
    """plt.figure()
    plt.scatter(x,y,c='r',marker='o') #numerical density
    plt.plot(z,g) #analytical gaussian function
    plt.xlabel("Rb/a")
    plt.ylabel("Density")
    plt.xlim(0.,2.0)
    plt.title("Density as a function of Rb (dimensionless)")
    plt.savefig("gaussian_%s.png" % i)"""


    #OPTICAL DEPTH new dataset
    a = []
    b = []

    dataset = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/optical_depth_NR=%s.txt' % i, 'r')

    for line in dataset:
        line = line.strip()
        A,B = line.split(',')
        a.append(float(A))
        b.append(float(B))

    c = [j-mean for j in a] #a-mean
    d = np.asarray((1.75668e+15/2.)*erf(c/(np.sqrt(2)*sd))+(1.75668e+15/2.)) #analytical error function
    ba = np.asarray(b) #b as an array
    e = (ba-d)/d #fractional error of numerical and analytical solutions
    rms = np.sqrt(e**2) #rms of fractional error
    err = average(rms) #average rms for each NR

    error.append(float(err))

    #plotting optical depth as a function of Rb
    dataset.close()
    """plt.figure()
    plt.scatter(a,b,c='r',marker='o') #numerical optical depth
    plt.plot(z,f) #analytical erf function
    plt.xlabel("Rb/a")
    plt.ylabel("Optical Depth")
    plt.title("Optical depth as a function of Rb (dimensionless)")
    plt.savefig("optical_depth_%s.png" % i)"""

#creating NR vs error textfile
dat = np.array([NR, error])
dat=dat.T
np.savetxt("/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/NR_error.txt", dat, fmt="%.0f,%.10f", delimiter = ',')

#plotting error(NR)
plt.figure()
plt.scatter(NR,error,c='r', marker='o')
plt.xlabel("NR")
plt.ylabel("error")
plt.savefig("error.png")
