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

fig = plt.figure(0,figsize=(10,4))
# fig.subplots_adjust(left=0.1, right =0.95,\
#         bottom=0.15, top =0.95, wspace=0.25)
sp1 = fig.add_subplot(121)
sp2 = fig.add_subplot(122)


# T = 5000     # temperature
dataDir = 'data/rg/lj/'


N = 1253     # bead number
# data from simulation
for T in [1, 2, 5, 10, 20, 50, \
        100, 200, 500, 1000, 2000, 5000]: 
    fname = dataDir + 'r_N'+str(N)+'_T'\
            +str(T)+'_5489.dat'
    data = np.loadtxt(fname, comments = '#')
    beadPos = data.reshape([-1,N,3])
    S = GyrationTensorNumerical(beadPos)
    kappa, sigma = Asphericity(S)
    sp1.scatter(T, kappa,\
            marker = 'o', color = 'g')
    sp2.scatter(T, sigma,\
            marker = 'o', color = 'g')

# # theory
N = 280
Tlog = np.logspace(0, 5, 50)
kappa_t = np.zeros(len(Tlog))
sigma_t = np.zeros(len(Tlog))
count = 0
for T in Tlog:
    S = GyrationTensor(N, T)
    kappa_t[count], sigma_t[count] = Asphericity(S)
    count = count + 1

sp1.plot(Tlog, kappa_t, 'b-')
sp2.plot(Tlog, sigma_t, 'b-')

sp1.set_xscale('log')
sp1.set_xlim([1,1e+4])
sp1.set_ylim([0,1])
sp1.set_xlabel(r'$T$')
sp1.set_ylabel(r'$\left<\Delta\right>$')

sp2.set_xscale('log')
sp2.set_xlim([1,1e+4])
sp2.set_ylim([-1.2,1.2])
sp2.set_xlabel(r'$T$')
sp2.set_ylabel(r'$\left<\Sigma\right>$')

fig.tight_layout()
plt.show()
