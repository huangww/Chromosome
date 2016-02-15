import numpy as np
import matplotlib.pyplot as plt
from findTau1D import *

font = {'family' : 'scans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 12}
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text', usetex = True)

figDir = 'fig/'
fig = plt.figure(0, figsize=(12,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

colors = ['b','r','g']
Temp = np.array([5, 10, 20])
for T,c in zip(Temp, colors):
    t, watchPara = GetWatchPara(T)
    # t, watchPara = LoadWatchPara(T)
    ax1.plot(t, watchPara,'.', color=c)
    # ax.plot(t, watchPara,'.', color=c)

    x, y = FindLinear(t, np.log(watchPara))
    coefs = np.polyfit(x,y,1)
    tau = -1.0/coefs[0]
    yFit = t * coefs[0] + coefs[1]
    tau = -1.0/coefs[0]
    labelString = 'T=' + str(T).zfill(2) + ', tau=' + str(int(tau))
    ax1.plot(t, np.exp(yFit), color=c, label=labelString)

Temp = np.array([1,2,5,10,20,50,100])
for T in Temp:
    t, watchPara = GetWatchPara(T)
    # t, watchPara = LoadWatchPara(T)
    x, y = FindLinear(t, np.log(watchPara))
    coefs = np.polyfit(x,y,1)
    tau = -1.0/coefs[0]
    ax2.scatter(T, tau)

ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r"$d_z$")
ax1.set_yscale('log')
ax1.legend(loc='lower left')
ax2.set_xlabel('$\\tilde{T}$')
ax2.set_ylabel(r"$\tau$")
plt.savefig(figDir + 'TempTau.pdf')
plt.show()
