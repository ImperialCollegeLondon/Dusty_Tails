import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import math
import pprint
from scipy.special import erf, erfc
from scipy import optimize
from textwrap import wrap

plt.ioff()

a = np.zeros(24)

NT = 4
NP = 3
NR = 2

a1 = []
a2 = []
a3 = []
den = [] * (NT*NP*NR)

b = np.reshape(a,(NT,NP,NR),order='C')

print(np.shape(b))
