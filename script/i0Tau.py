import numpy as np
import matplotlib.pyplot as plt
# from findTau import *
from findTau1D import *

font = {'family' : 'scans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 16}
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text', usetex = True)
plt.close('all')
fig = plt.figure(0, figsize=(12,5))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
fig.subplots_adjust(left=0.12, right=0.95,
        bottom=0.15, top=0.95, wspace=0.25)

N = 100
T = 1

figDir = 'fig/'
dataDir = 'data/'
fname = dataDir+'asep_mean_N'+str(N)+'_T'+str(T)+'.dat'
data = np.loadtxt(fname)
# fname = dataDir+'mean_N'+str(N)+'_T'+str(T)+'.dat'
# data = np.loadtxt(fname).reshape([-1,N,3])

toReal = 18.5

i0 = [1,10, 20, 50]
# colors = ['b','r','y','m','g']
for i in i0:
    t, watchPara = GetWatchPara(data, i0=i)
    watchPara = watchPara/watchPara.max()
    ax1.plot(t*toReal, watchPara,'-',label=r'$i=$'+str(i))
    markers, = ax2.plot(t*toReal, watchPara,'-')

    # x, y = FindLinear(i, t, np.log(np.maximum(watchPara, 1e-9, watchPara)))
    # # x, y = t[t0:tend], np.log(np.maximum(watchPara[t0:tend], 1e-9, watchPara[t0:tend]))
    # coefs = np.polyfit(x,y,1)
    # yFit = t * coefs[0] + coefs[1]
    # tau = -1.0/coefs[0]
    # labelString = r'$i=$ ' + str(i).zfill(2) + ', '+ r'$\tau=$ ' + str(int(tau*toReal))
    # ax2.plot(t*toReal, np.exp(yFit), color=markers.get_color(), label=labelString)

ax1.set_xlabel(r'$t[s]$')
ax1.set_ylabel(r"$d_z$")
ax1.legend()
ax2.set_xlabel(r'$t[s]$')
ax2.set_ylabel(r"$d_z$")
ax2.set_yscale('log')
ax1.set_xlim([0,15000])
ax2.set_xlim([0,15000])
# ax2.set_ylim(1e-25)
# ax2.legend(loc='lower left')
fname = 'I0Tau_N'+str(N)+'_T' + str(T) + '.pdf'
fig.savefig(figDir + fname)
plt.show()
