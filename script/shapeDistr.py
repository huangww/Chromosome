import numpy as np
import matplotlib.pyplot as plt
from theoryGyration import *
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

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

N = 1253
dataDir = 'data/rg/lj/'

T = 1
fname = dataDir + 'r_N'+str(N)+'_T'\
        +str(T)+'_5489.dat'
data = np.loadtxt(fname, comments = '#')
beadPos = data.reshape([-1,N,3])
kappa, sigma = AsphericityArray(beadPos)
rho = 2*np.sqrt(kappa)
theta = np.arccos(sigma)/3
sp1.hist2d(theta, rho, bins=50, normed=True,\
        range=np.array([(0,np.pi/3),(0,2)]),\
        cmap='YlOrRd')
sp1.set_xlim([0,np.pi/3])
sp1.set_xticks([0,np.pi/6,np.pi/3])
sp1.set_xticklabels([0,r"$\pi/6$",r"$\pi/3$"])
sp1.set_ylim([0,2])
sp1.set_xlabel(r'$\theta$')
sp1.set_ylabel(r"$\rho$")
# sp1.legend([r"T="+str(T)])
axins = zoomed_inset_axes(sp1, 30, loc=10) # zoom = 6
axins.hist2d(theta, rho, bins=50, normed=True,\
        range=np.array([(0,0.02),(1.95,2)]), \
        cmap='YlOrRd')
# axins.set_xlim(0.5, 1.0)
# axins.set_ylim(0.5, 1.0)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(sp1, axins, loc1=3, loc2=1, fc="none", ec="0.5")


T = 10
fname = dataDir + 'r_N'+str(N)+'_T'\
        +str(T)+'_5489.dat'
data = np.loadtxt(fname, comments = '#')
beadPos = data.reshape([-1,N,3])
kappa, sigma = AsphericityArray(beadPos)
rho = 2*np.sqrt(kappa)
theta = np.arccos(sigma)/3
sp2.hist2d(theta, rho, bins=50, normed=True,\
        range=np.array([(0,np.pi/3),(0,2)]),\
        cmap='YlOrRd')
sp2.set_xlim([0,np.pi/3])
sp2.set_xticks([0,np.pi/6,np.pi/3])
sp2.set_xticklabels([0,r"$\pi/6$",r"$\pi/3$"])
sp2.set_ylim([0,2])
sp2.set_xlabel(r'$\theta$')
sp2.set_ylabel(r"$\rho$")
# sp2.legend([r"T="+str(T)])
axins = zoomed_inset_axes(sp2, 6, loc=10) # zoom = 6
axins.hist2d(theta, rho, bins=50, normed=True,\
        range=np.array([(0,0.1),(1.75,2)]), \
        cmap='YlOrRd')
# axins.set_xlim(0.5, 1.0)
# axins.set_ylim(0.5, 1.0)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(sp2, axins, loc1=3, loc2=1, fc="none", ec="0.5")


T = 50
fname = dataDir + 'r_N'+str(N)+'_T'\
        +str(T)+'_5489.dat'
data = np.loadtxt(fname, comments = '#')
beadPos = data.reshape([-1,N,3])
kappa, sigma = AsphericityArray(beadPos)
rho = 2*np.sqrt(kappa)
theta = np.arccos(sigma)/3
sp3.hist2d(theta, rho, bins=50, normed=True,\
        range=np.array([(0,np.pi/3),(0,2)]),\
        cmap='YlOrRd')
sp3.set_xlim([0,np.pi/3])
sp3.set_xticks([0,np.pi/6,np.pi/3])
sp3.set_xticklabels([0,r"$\pi/6$",r"$\pi/3$"])
sp3.set_ylim([0,2])
sp3.set_xlabel(r'$\theta$')
sp3.set_ylabel(r"$\rho$")
# sp3.legend([r"T="+str(T)])

T = 5000
fname = dataDir + 'r_N'+str(N)+'_T'\
        +str(T)+'_5489.dat'
data = np.loadtxt(fname, comments = '#')
beadPos = data.reshape([-1,N,3])
kappa, sigma = AsphericityArray(beadPos)
rho = 2*np.sqrt(kappa)
theta = np.arccos(sigma)/3
sp4.hist2d(theta, rho, bins=50, normed=True,\
       range=np.array([(0,np.pi/3),(0,2)]),\
        cmap='YlOrRd')
sp4.set_xlim([0,np.pi/3])
sp4.set_xticks([0,np.pi/6,np.pi/3])
sp4.set_xticklabels([0,r"$\pi/6$",r"$\pi/3$"])
sp4.set_ylim([0,2])
sp4.set_xlabel(r'$\theta$')
sp4.set_ylabel(r"$\rho$")
# sp4.legend([r"T="+str(T)])

fig.tight_layout()
plt.show()
