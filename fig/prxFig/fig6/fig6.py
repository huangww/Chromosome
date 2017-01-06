import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

def tauTheory(a, b, L):
    tau = 1. / (a + b - 2.*np.sqrt(a*b)*np.cos(np.pi/L))
    return tau

fig = plt.figure(0,figsize=(5, 3.6))
plt.rc('text', usetex=True)
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 15 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
fig.subplots_adjust(left=0.13, right =0.95,\
        bottom=0.13, top =0.95, wspace=0.25)

ax = fig.add_subplot(111)
L = 100
tau0 = (L/np.pi)**2/3
Teff = np.linspace(1, 1000, 1000)
factor = np.exp(-1./Teff)
a = 6.*factor/(1.+factor)
b = 6./(1.+factor)
ax.plot(Teff, tauTheory(a, b, L)/tau0, 'b-')
fname = 'BeadRod_tau_N100.dat'
tauData = np.loadtxt(fname)
ax.plot(tauData[:,0], tauData[:,-1]/tau0, 'bo')
# mark the data point will show acf in  inset

ax.set_xscale('log')
ax.set_xlabel(r'$k_B T/\Delta E$')
ax.set_ylabel(r'$\tau_{\parallel}/\tau_{Rouse}$')
ax.legend(['Theory', 'Simulation'], loc='upper left', fontsize=15)
ax.set_xlim([0,1000])
ax.set_ylim([0,1.2])

axinset = plt.axes([0.6, 0.26, 0.32, 0.3])
axinset.plot(Teff, tauTheory(a, b, L)/tau0, 'b-')
fname = 'BeadRod_tau_N100_p.dat'
tauData = np.loadtxt(fname)
axinset.plot(tauData[:,0], tauData[:,-1]/tau0, 'bo')
axinset.set_xlabel(r'$k_B T/\Delta E$', labelpad=0)
axinset.set_ylabel(r'$\tau_{\perp} / \tau_{Rouse}$',labelpad=0)
axinset.set_xscale('log')
axinset.set_xlim([0,100])
axinset.set_ylim([0,1.2])
axinset.set_yticks([0, 0.4, 0.8, 1.2])
axinset.tick_params(axis='both', which='major', pad=2)

plt.show()
