import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
from scipy.special import erf, erfc
from scipy import optimize
from textwrap import wrap

plt.ioff()

#defining variables
rho = 1.4654e-14 #max density in density gaussian (dimensionless)
kappa = 7.24924e+13 #opacity as calculated in cpp file (dimensionless)
mean = 1.
variance = (0.2*mean)**2
sd = math.sqrt(variance)
z = np.linspace(0., 2.0, 1000.)
g = rho*(np.sqrt(2.*math.pi)*sd)*stats.norm.pdf(z,mean,sd) #analytic gaussian curve
f = (np.sqrt(math.pi/2.)*sd*(rho*kappa))*(erf((z-mean)/(np.sqrt(2)*sd))+1) #analytic erf function +((np.sqrt(math.pi/2.)*sd*(1.4654e-14*7.24924e+13))/2.)

#defining functions
def average(x): #finds the average of a list
    return sum(x)/len(x)

def factor_two(n):
    while n < 1000:
        yield n
        n += 1

#defining lists to plot error against NR
error = []
logerror = []
NR = []
logNR = []

#for loop to produce optical depth plots for different NR
for i in factor_two(10):

    #DENSITY construct dataset to plot
    x = []
    y = []
    NR.append(i)
    logNR.append(math.log(float(i)))

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
    plt.ylim(0.,3.e-14)
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

    #plotting optical depth as a function of Rb
    """dataset.close()
    plt.figure()
    plt.scatter(a,b,c='r',marker='o') #numerical optical depth
    plt.plot(z,f) #analytical erf function
    plt.xlabel("Rb/a")
    plt.ylabel("Optical Depth")
    plt.title("Optical depth as a function of Rb (dimensionless)")
    plt.savefig("optical_depth_%s.png" % i)"""

    #error calculations
    c = [j-mean for j in a] #a-mean
    d = np.asarray((np.sqrt(math.pi/2.)*sd*(rho*kappa))*(erf((c)/(np.sqrt(2)*sd))+1)) #analytical error function
    ba = np.asarray(b) #b as an array
    e = (ba-d)/d #fractional error of numerical and analytical solutions
    rms = np.sqrt(e**2) #rms of fractional error
    err = average(rms) #average rms for each NR

    error.append(float(err))
    logerror.append(math.log(float(err)))


#creating NR vs error textfile
dat = np.array([NR, error])
dat=dat.T
np.savetxt("/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/NR_error.txt", dat, fmt="%.0f,%.10f", delimiter = ',')

#plotting error(NR)
plt.figure()
plt.scatter(NR,error,c='r', marker='o')
plt.xlabel("NR")
plt.ylabel("Error")
plt.title("\n".join(wrap("Error between the analytical erf and numerical solution for optical depth as a function of NR, the number of grid cells")))
plt.savefig("error.png")

plt.figure()
plt.scatter(logNR,logerror,c='r',marker='o')
plt.xlabel("log(NR)")
plt.ylabel("log(error)")
plt.title("\n".join(wrap("A log plot of the error between the analytical erf and numerical solution for optical depth as a function of NR, the number of grid cells")))
plt.savefig("error_logplot.png")
