import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

def tauTheory(a, b, L):
    tau = 1. / (a + b - 2.*np.sqrt(a*b)*np.cos(np.pi/L))
    return tau

fig = plt.figure(0,figsize=(6.68, 4.46))
plt.rc('text', usetex=True)
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 15 }
plt.rc('lines', lw=2)
plt.rc('font', **font)

ax = fig.add_subplot(111)
L = 100
tau0 = (L/np.pi)**2
Teff = np.linspace(1, 1000, 1000)
factor = np.exp(-1./Teff)
a = 2.*factor/(1.+factor)
b = 2./(1.+factor)
ax.plot(Teff, tauTheory(a, b, L)/tau0, '-')
fname = 'Asep_tau_N100_0.dat'
tauData = np.loadtxt(fname)
ax.plot(tauData[:,0], tauData[:,-1]/tau0, 'o')
# mark the data point will show acf in  inset

ax.set_xscale('log')
ax.set_xlabel(r'$\log(\beta/\alpha)$')
ax.set_ylabel(r'$\tau / \tau_{Rouse}$')
ax.legend(['Theory', 'Simulation'], loc='upper left')
ax.set_xlim([0,1000])
ax.set_ylim([0,1.2])

plt.show()
