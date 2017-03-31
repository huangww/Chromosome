import numpy as np
import matplotlib.pyplot as plt
import acfFit as af
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

def tauTheory(a, b, L):
    tau = 1. / (a + b - 2.*np.sqrt(a*b)*np.cos(np.pi/L))
    return tau

fig = plt.figure(0,figsize=(5, 3.6))
plt.rc('text', usetex=True)
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 14 }
plt.rc('lines', lw=2)
plt.rc('font', **font)

fig.subplots_adjust(left=0.13, right =0.95,\
        bottom=0.13, top =0.95, wspace=0.25)

ax = fig.add_subplot(111)
L = 100
tau0 = (L/np.pi)**2
Teff = np.linspace(1, 1000, 1000)
factor = np.exp(-1./Teff)
a = 2.*factor/(1.+factor)
b = 2./(1.+factor)
ax.plot(np.log(b/a), tauTheory(a, b, L)/tau0, 'b-')
fname = 'Asep_tau_N100.dat'
tauData = np.loadtxt(fname)
ax.plot(1./tauData[:,0], tauData[:,-1]/tau0, 'bo')
# mark the data point will show acf in  inset
for x in tauData:
    if x[0]==10:
        ax.plot(1./x[0], x[-1]/tau0, marker='o',\
                markersize=9, color='r', alpha=0.5)
    if x[0]==20:
        ax.plot(1./x[0], x[-1]/tau0, marker='o',\
                markersize=9, color='c', alpha=0.5)
    if x[0]==50:
        ax.plot(1./x[0], x[-1]/tau0, marker='o',\
                markersize=9, color='m', alpha=0.5)

ax.set_xscale('log')
ax.set_xlabel(r'$q$')
ax.set_ylabel(r'$\tau / \tau_{Rouse}$')
ax.legend(['Theory', 'Simulation'], loc='upper right', fontsize=15)
ax.set_xlim([0,1])
ax.set_ylim([0,1.2])


# axinset = inset_axes(ax, width="40%", height="40%", loc=7)
axinset = plt.axes([0.25, 0.26, 0.32, 0.3])
Tset = [10, 20, 50]
colors = ['r', 'c', 'm']
for T,c in zip(Tset,colors):
    fname = 'Asep_acf_N100_T%g.dat'%T
    t, acf = np.loadtxt(fname)
    axinset.plot(t, acf, '-', color=c)
    start, end = af.FindFitRange(acf)
    coefs = np.polyfit(t[start:end], np.log(acf[start:end]), 1)
    axinset.plot(t, np.exp(coefs[0]*t+coefs[1]), 'k--',lw=1)

axinset.set_yscale('log')
axinset.set_xlabel(r'$\Delta t$', labelpad=0)
axinset.set_ylabel('ACF',labelpad=0)
axinset.set_xlim([0, 4e3])
axinset.set_xticks([0, 1e3, 2e3, 3e3, 4e3])
axinset.ticklabel_format(axis='x', style='sci',scilimits=(0,1))
axinset.tick_params(axis='both', which='major', pad=1)
axinset.set_ylim([1e-2, 1])

plt.show()
