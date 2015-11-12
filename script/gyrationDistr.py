import numpy as np
import matplotlib.pyplot as plt
from theoryGyration import *

font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 16 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text',usetex = True)

fig = plt.figure(0,figsize=(12,9))
sp1 = fig.add_subplot(221)
sp2 = fig.add_subplot(222)
sp3 = fig.add_subplot(223)
sp4 = fig.add_subplot(224)

beadNumber = 100
dataDir = 'data/rg/'

T = 1
fileName = dataDir + 'rg_N' + str(beadNumber) + \
        '_T'+str(T) + '_5489.dat'
data = np.loadtxt(fileName,comments='#')
sp1.hist(data,bins=50,normed=True)
sp1.set_xlabel(r'$R_g^2$')
sp1.set_ylabel(r"$P(R_g^2)$")
sp1.legend(['T = ' + str(T)])

T = 10
fileName = dataDir + 'rg_N' + str(beadNumber) + \
        '_T'+str(T) + '_5489.dat'
data = np.loadtxt(fileName,comments='#')
sp2.hist(data,bins=50,normed=True)
sp2.set_xlabel(r'$R_g^2$')
sp2.set_ylabel(r"$P(R_g^2)$")
sp2.legend(['T = ' + str(T)])

T = 50
fileName = dataDir + 'rg_N' + str(beadNumber) + \
        '_T'+str(T) + '_5489.dat'
data = np.loadtxt(fileName,comments='#')
sp3.hist(data,bins=50,normed=True)
sp3.set_xlabel(r'$R_g^2$')
sp3.set_ylabel(r"$P(R_g^2)$")
sp3.legend(['T = ' + str(T)])

T = 5000
meanRgs = beadNumber/12.0
t = np.array([0.1, 0.2, 0.3, 0.4, 0.5, \
    0.6, 0.7, 0.8, 0.9, 1.0, \
    1.2, 1.4, 1.6, 1.8, 2.0, \
    2.4,2.8])
prgs = np.array([0.0, 0.0, 0.0143, 0.1535, 0.4878, \
    0.8733, 1.1462, 1.2498, 1.2100, 1.0804, \
    0.7345, 0.4359, 0.2373, 0.1215, 0.0596, \
    0.0131, 0.0026])
fileName = dataDir + 'rg_N' + str(beadNumber) + \
        '_T'+str(T) + '_5489.dat'
data = np.loadtxt(fileName,comments='#')
sp4.hist(data,bins=50,normed=True)
sp4.plot(t*meanRgs, prgs/meanRgs, 'r-')
sp4.set_xlabel(r'$R_g^2$')
sp4.set_ylabel(r"$P(R_g^2)$")
sp4.legend(['High $T$ asymptotic','T = ' + str(T)])

fig.tight_layout()
plt.show()
