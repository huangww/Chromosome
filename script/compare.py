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

toReal = 1

i = 50
fname = dataDir+'asep_mean_N'+str(N)+'_T'+str(T)+'.dat'
data = np.loadtxt(fname)
t, watchPara = GetWatchPara(data, i0=i)
watchPara = watchPara/watchPara.max()
ax1.plot(t, watchPara,'.-')
ax2.plot(t, watchPara,'.-')

fname = dataDir+'asep_var_N'+str(N)+'_T'+str(T)+'.dat'
data = np.loadtxt(fname)
t, watchPara = GetWatchPara(data, i0=i)
watchPara = watchPara/watchPara.max()
ax1.plot(t, watchPara,'.-')
ax2.plot(t, watchPara,'.-')

ax1.set_xlabel(r'$t$')
ax1.set_ylabel(r"$d_z$")
ax2.set_xlabel(r'$t$')
ax2.set_ylabel(r"$d_z$")
ax2.set_yscale('log')
ax1.set_xlim([0,800])
ax2.set_xlim([0,800])
fname = 'compairMeanVar_N'+str(N)+'_T' + str(T) + '.pdf'
fig.savefig(figDir + fname)
plt.show()
