import numpy as np
import matplotlib.pyplot as plt
import statsmodels.tsa.stattools as ss


plt.close('all')
fig = plt.figure(0)


ax = fig.add_subplot(131)
fname = 'data/rg1D_N100_T1_eq.dat'
data = np.loadtxt(fname)
nlags = next(i for i,v in enumerate(data[:,0]) if v>1e2)
acfRg = ss.acf(data[:,1], nlags=nlags, fft=True)
acfR1 = ss.acf(data[:,2], nlags=nlags, fft=True)
acfRd = ss.acf(data[:,3], nlags=nlags, fft=True)
t = data[:nlags+1,0]
ax.plot(t, acfRg, label='$R_g$')
ax.plot(t, acfR1, label='$R_1$')
ax.plot(t, acfRd, label='$R_d$')
ax.legend()
ax.set_xlabel(r'$\delta t$')
ax.set_ylabel(r'$ACF$')
ax.set_xlim(0,1e2)
ax.set_yscale('log')
ax.set_title(r'$1/F=1$')

ax = fig.add_subplot(132)
fname = 'data/rg1D_N100_T100_eq.dat'
data = np.loadtxt(fname)
nlags = next(i for i,v in enumerate(data[:,0]) if v>1e4)
acfRg = ss.acf(data[:,1], nlags=nlags, fft=True)
acfR1 = ss.acf(data[:,2], nlags=nlags, fft=True)
acfRd = ss.acf(data[:,3], nlags=nlags, fft=True)
t = data[:nlags+1,0]
ax.plot(t, acfRg, label='$R_g$')
ax.plot(t, acfR1, label='$R_1$')
ax.plot(t, acfRd, label='$R_d$')
# ax.legend()
ax.set_xlabel(r'$\delta t$')
ax.set_ylabel(r'$ACF$')
ax.set_xlim(0,1e4)
ax.set_yscale('log')
ax.set_title(r'$1/F=100$')

ax = fig.add_subplot(133)
fname = 'data/rg1D_N100_T1000_eq.dat'
data = np.loadtxt(fname)
nlags = next(i for i,v in enumerate(data[:,0]) if v>1e4)
acfRg = ss.acf(data[:,1], nlags=nlags, fft=True)
acfR1 = ss.acf(data[:,2], nlags=nlags, fft=True)
acfRd = ss.acf(data[:,3], nlags=nlags, fft=True)
t = data[:nlags+1,0]
ax.plot(t, acfRg, label='$R_g$')
ax.plot(t, acfR1, label='$R_1$')
ax.plot(t, acfRd, label='$R_d$')
# ax.legend()
ax.set_xlabel(r'$\delta t$')
ax.set_ylabel(r'$ACF$')
ax.set_xlim(0,1e4)
ax.set_yscale('log')
ax.set_title(r'$1/F=1000$')

fname = 'fig/checkLongest.pdf'
# fig.savefig(fname)
plt.show()

