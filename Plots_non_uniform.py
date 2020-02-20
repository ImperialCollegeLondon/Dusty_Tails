import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
from scipy.special import erf, erfc
from scipy import optimize
from textwrap import wrap

plt.ioff()

#defining variables
mean = 1.
variance = (0.2*mean)**2
sd = math.sqrt(variance)
z = np.linspace(0., 2.0, 1000.)
g = 1.469688471e-14*(np.sqrt(2*np.pi)*sd)*stats.norm.pdf(z,mean,sd) #analytic gaussian curve
f = ((1.0715748/2.)*(sd*np.sqrt(np.pi*2.)))*(erf((z-mean)/(np.sqrt(2)*sd))+1.) #analytic erf function

#defining functions
def average(x): #finds the average of a list
    return sum(x)/len(x)

def factor_two(n):
    while n < 101:
        yield n
        n += 1

#defining lists to plot error against NR
error = []
erroru = []
NR = []
NRu = []
logerror = []
logerroru = []
logNR = []
logNRu = []

#UNIFORM GRID

#for loop to produce optical depth plots for different NR
for i in factor_two(100):

    #DENSITY (non-uniform grid) construct dataset to plot
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
    plt.figure()
    plt.scatter(x,y,c='r',marker='o') #numerical density
    plt.plot(z,g) #analytical gaussian function
    plt.xlabel("Rb/a")
    plt.ylabel("Density")
    plt.xlim(0.,2.0)
    plt.ylim(0.,2.0e-14)
    plt.title("Density as a function of Rb (uniform grid) for NR=100")
    plt.savefig("gaussianu_%s.png" % i)


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
    dataset.close()
    plt.figure()
    plt.scatter(a,b,c='r',marker='o') #numerical optical depth
    plt.plot(z,f) #analytical erf function
    plt.xlabel("Ra/a")
    plt.ylabel("Optical Depth")
    plt.title("Optical depth as a function of Ra (uniform grid) for NR=100")
    plt.savefig("optical_depthu_%s.png" % i)

    c = [j-mean for j in a] #a-mean
    d = np.asarray((1.0715748/2.)*(sd*np.sqrt(np.pi*2.))*(erf(c/(np.sqrt(2)*sd))+1.)) #analytical error function
    ba = np.asarray(b) #b as an array
    e = (ba-d)/d #fractional error of numerical and analytical solutions
    rms = np.sqrt(e**2) #rms of fractional error
    erru = average(rms) #average rms for each NR

    erroru.append(float(erru))
    logerroru.append(math.log(float(erru)))

#NON UNIFORM GRID

    #DENSITY (non-uniform grid) construct dataset to plot
    p = []
    q = []
    NR.append(i)
    logNR.append(math.log(float(i)))

    dataset = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/gaussian_NR=%s.txt' % i,'r')

    for line in dataset:
        line = line.strip()
        P,Q = line.split(',')
        p.append(float(P))
        q.append(float(Q))

    dataset.close()
    plt.figure()
    plt.scatter(p,q,c='r',marker='o') #numerical density
    plt.plot(z,g) #analytical gaussian function
    plt.xlabel("Rb/a")
    plt.ylabel("Density")
    plt.xlim(0.,2.0)
    plt.ylim(0.,2.0e-14)
    plt.title("Density as a function of Rb (gaussian grid) for NR=100")
    plt.savefig("gaussian_%s.png" % i)


    #OPTICAL DEPTH new dataset
    r = []
    s = []

    dataset = open('/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/optical_depth_NR=%s.txt' % i, 'r')

    for line in dataset:
        line = line.strip()
        R,S = line.split(',')
        r.append(float(R))
        s.append(float(S))

    #plotting optical depth as a function of Rb
    dataset.close()
    plt.figure()
    plt.scatter(r,s,c='r',marker='o') #numerical optical depth
    plt.plot(z,f) #analytical erf function
    plt.xlabel("Ra/a")
    plt.ylabel("Optical Depth")
    plt.title("Optical depth as a function of Ra (gaussian grid) for NR=100")
    plt.savefig("optical_depth_%s.png" % i)

    C = [j-mean for j in r] #r-mean
    D = np.asarray((1.0715748/2.)*(sd*np.sqrt(np.pi*2.))*(erf(C/(np.sqrt(2)*sd))+1.)) #analytical error function
    sa = np.asarray(s) #s as an array
    E = (sa-D)/D #fractional error of numerical and analytical solutions
    rms = np.sqrt(E**2) #rms of fractional error
    err = average(rms) #average rms for each NR

    error.append(float(err))
    logerror.append(math.log(float(err)))

#creating NR vs error textfile uniform grid
datu = np.array([NRu, erroru])
datu = datu.T
np.savetxt("/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/NR_erroru.txt", datu, fmt="%.0f,%.10f", delimiter = ',')

#creating NR vs error textfile
dat = np.array([NR, error])
dat = dat.T
np.savetxt("/Users/annawilson/Documents/GitHub/Dusty_Tails/Text_files/NR_error.txt", dat, fmt="%.0f,%.10f", delimiter = ',')

#plotting error(NR) - uniform
plt.figure()
plt.scatter(NRu,erroru,c='b', marker='o')
plt.xlabel("NR")
plt.ylabel("error")
plt.title("\n".join(wrap("Error between the analytical erf and numerical solution for optical depth as a function of NR (uniform grid)")))
plt.savefig("error_uniform.png")

plt.figure()
plt.scatter(logNRu,logerroru,c='b',marker='o')
plt.xlabel("log(NR)")
plt.ylabel("log(error)")
plt.title("\n".join(wrap("A log plot of the error between the analytical erf and numerical solution for optical depth as a function of NR (uniform grid)")))
plt.savefig("error_logplot_uniform.png")

#plotting error(NR) - non-uniform
plt.figure()
plt.scatter(NR,error,c='r', marker='o')
plt.xlabel("NR")
plt.ylabel("error")
plt.title("\n".join(wrap("Error between the analytical erf and numerical solution for optical depth as a function of NR (gaussian grid)")))
plt.savefig("error.png")

plt.figure()
plt.scatter(logNR,logerror,c='r',marker='o')
plt.xlabel("log(NR)")
plt.ylabel("log(error)")
plt.title("\n".join(wrap("A log plot of the error between the analytical erf and numerical solution for optical depth as a function of NR (gaussian grid)")))
plt.savefig("error_logplot.png")

#plotting comparison of error for non-uniform and uniform grids
plt.figure()
plt.scatter(NR,error,c='r',marker='o')
plt.scatter(NRu,erroru,c='b',marker='o')
plt.xlabel("log(NR)")
plt.ylabel("log(error)")
plt.title("\n".join(wrap("Comparison of error between results from uniform (blue) and gaussian (red) grid")))
plt.savefig("error_comparison.png")

plt.figure()
plt.scatter(logNR,logerror,c='r',marker='o')
plt.scatter(logNRu,logerroru,c='b',marker='o')
plt.xlabel("log(NR)")
plt.ylabel("log(error)")
plt.title("\n".join(wrap("A log plot of the comparison of error between results from uniform (blue) and gaussian (red) grid")))
plt.savefig("error_log_comparison.png")
