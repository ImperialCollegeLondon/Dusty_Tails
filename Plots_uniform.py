import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
from scipy.special import erf, erfc
from scipy import optimize
from textwrap import wrap

plt.ioff()

#defining variables
"""rho = 1.4654e-14 #max density in density gaussian (dimensionless)
kappa = 7.24924e+13 #opacity as calculated in cpp file (dimensionless)
mean = 1.
variance = (0.2*mean)**2
sd = math.sqrt(variance)
z = np.linspace(0., 2.0, 1000.)
g = rho*(np.sqrt(2.*math.pi)*sd)*stats.norm.pdf(z,mean,sd) #analytic gaussian curve
f = (np.sqrt(math.pi/2.)*sd*(rho*kappa))*(erf((z-mean)/(np.sqrt(2)*sd))+1) #analytic erf function +((np.sqrt(math.pi/2.)*sd*(1.4654e-14*7.24924e+13))/2.)
"""

mean = 1.
variance = (0.2*mean)**2
sd = math.sqrt(variance)
z = np.linspace(0., 2.0, 1000.)
g = 1.44e-14*stats.norm.pdf(z,mean,sd) #analytic gaussian curve
f = (1.05/2.)*erf((z-mean)/(np.sqrt(2)*sd))+(1.05/2.) #analytic erf function
"""g = 40.3877*stats.norm.pdf(z,mean,sd) #analytic gaussian curve
f = (1.75668e+15/2.)*erf((z-mean)/(np.sqrt(2)*sd))+(1.75668e+15/2.) #analytic erf function"""

#defining functions
def average(x): #finds the average of a list
    return sum(x)/len(x)

def factor_two(n):
    while n < 1002:
        yield n
        n += 1

#defining lists to plot error against NR
erroru = []
logerroru = []
NRu = []
logNRu = []

#for loop to produce optical depth plots for different NR
for i in factor_two(11):

    #DENSITY construct dataset to plot
    x = []
    y = []
    NRu.append(i)
    logNRu.append(math.log(float(i)))

    dataset = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/gaussian_NRu=%s.txt' % i,'r')

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

    dataset = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/optical_depth_NRu=%s.txt' % i, 'r')

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
    d = np.asarray((1.75668e+15/2.)*erf(c/(np.sqrt(2)*sd))+(1.75668e+15/2.))
    """d = np.asarray((np.sqrt(math.pi/2.)*sd*(rho*kappa))*(erf((c)/(np.sqrt(2)*sd))+1))""" #analytical error function
    ba = np.asarray(b) #b as an array
    e = (ba-d)/d #fractional error of numerical and analytical solutions
    rms = np.sqrt(e**2) #rms of fractional error
    erru = average(rms) #average rms for each NR

    erroru.append(float(erru))
    logerroru.append(math.log(float(erru)))


#creating NR vs error textfile
datu = np.array([NRu, erroru])
datu=datu.T
np.savetxt("/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/NR_erroru.txt", datu, fmt="%.0f,%.10f", delimiter = ',')

#plotting error(NR)
plt.figure()
plt.scatter(NRu,erroru,c='r', marker='o')
plt.xlabel("NR")
plt.ylabel("Error")
plt.title("\n".join(wrap("Error between the analytical erf and numerical solution for optical depth as a function of NR, the number of grid cells")))
plt.savefig("error_uniform.png")

plt.figure()
plt.scatter(logNRu,logerroru,c='r',marker='o')
plt.xlabel("log(NR)")
plt.ylabel("log(error)")
plt.title("\n".join(wrap("A log plot of the error between the analytical erf and numerical solution for optical depth as a function of NR, the number of grid cells")))
plt.savefig("error_logplot_uniform.png")
