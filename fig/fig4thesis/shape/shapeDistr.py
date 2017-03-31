import numpy as np
import matplotlib.pyplot as plt
from theoryGyration import *
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 16 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text',usetex = True)

fig = plt.figure(0,figsize=(8,6))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
fig.subplots_adjust(left=0.1, right =0.95,\
        bottom=0.1, top =0.95, wspace=0.25)

N = 100

T = 1
fname = 'BeadRod_r_N'+str(N)+'_F%g_T1_0.dat'%(1./T)
data = np.loadtxt(fname, comments = '#')
beadPos = data.reshape([-1,N,3])
kappa, sigma = AsphericityArray(beadPos)
rho = 2*np.sqrt(kappa)
theta = np.arccos(sigma)/3
hax1 = ax1.hist2d(theta, rho, bins=50, normed=True,\
        range=np.array([(0,np.pi/3),(0,2)]),\
        cmap='YlOrRd')
ax1.set_xlim([0,np.pi/3])
ax1.set_xticks([0,np.pi/6,np.pi/3])
ax1.set_xticklabels([r'$0$',r'$\pi/6$',r'$\pi/3$'])
ax1.set_ylim([0,2])
ax1.set_xlabel(r'$\theta$',labelpad=0)
ax1.set_ylabel(r'$\rho$',labelpad=0)
axColorbar = inset_axes(ax1, width='3%',  height='30%', loc=1)
fig.colorbar(hax1[3], cax=axColorbar, ticks=[])
# ax1.legend([r'T='+str(T)])
axins = zoomed_inset_axes(ax1, 30, loc=10) # zoom = 6
axins.hist2d(theta, rho, bins=50, normed=True,\
        range=np.array([(0,0.02),(1.95,2)]), \
        cmap='YlOrRd')
# axins.set_xlim(0.5, 1.0)
# axins.set_ylim(0.5, 1.0)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax1, axins, loc1=3, loc2=1, fc='none', ec='0.5')


T = 10
fname = 'BeadRod_r_N'+str(N)+'_F%g_T1_0.dat'%(1./T)
data = np.loadtxt(fname, comments = '#')
beadPos = data.reshape([-1,N,3])
kappa, sigma = AsphericityArray(beadPos)
rho = 2*np.sqrt(kappa)
theta = np.arccos(sigma)/3
hax2 = ax2.hist2d(theta, rho, bins=50, normed=True,\
        range=np.array([(0,np.pi/3),(0,2)]),\
        cmap='YlOrRd')
ax2.set_xlim([0,np.pi/3])
ax2.set_xticks([0,np.pi/6,np.pi/3])
ax2.set_xticklabels([r'$0$',r'$\pi/6$',r'$\pi/3$'])
ax2.set_ylim([0,2])
ax2.set_xlabel(r'$\theta$',labelpad=0)
ax2.set_ylabel(r'$\rho$',labelpad=0)
axColorbar = inset_axes(ax2, width='3%',  height='30%', loc=1)
fig.colorbar(hax2[3], cax=axColorbar, ticks=[])
# ax2.legend([r'T='+str(T)])
axins = zoomed_inset_axes(ax2, 6, loc=10) # zoom = 6
axins.hist2d(theta, rho, bins=50, normed=True,\
        range=np.array([(0,0.1),(1.75,2)]), \
        cmap='YlOrRd')
# axins.set_xlim(0.5, 1.0)
# axins.set_ylim(0.5, 1.0)
plt.xticks(visible=False)
plt.yticks(visible=False)
mark_inset(ax2, axins, loc1=3, loc2=1, fc='none', ec='0.5')


T = 100
fname = 'BeadRod_r_N'+str(N)+'_F%g_T1_0.dat'%(1./T)
data = np.loadtxt(fname, comments = '#')
beadPos = data.reshape([-1,N,3])
kappa, sigma = AsphericityArray(beadPos)
rho = 2*np.sqrt(kappa)
theta = np.arccos(sigma)/3
hax3 = ax3.hist2d(theta, rho, bins=50, normed=True,\
        range=np.array([(0,np.pi/3),(0,2)]),\
        cmap='YlOrRd')
ax3.set_xlim([0,np.pi/3])
ax3.set_xticks([0,np.pi/6,np.pi/3])
ax3.set_xticklabels([r'$0$',r'$\pi/6$',r'$\pi/3$'])
ax3.set_ylim([0,2])
ax3.set_xlabel(r'$\theta$',labelpad=0)
ax3.set_ylabel(r'$\rho$',labelpad=0)
# ax3.legend([r'T='+str(T)])
axColorbar = inset_axes(ax3, width='3%',  height='30%', loc=1)
fig.colorbar(hax3[3], cax=axColorbar, ticks=[])

T = 1000
fname = 'BeadRod_r_N'+str(N)+'_F%g_T1_0.dat'%(1./T)
data = np.loadtxt(fname, comments = '#')
beadPos = data.reshape([-1,N,3])
kappa, sigma = AsphericityArray(beadPos)
rho = 2*np.sqrt(kappa)
theta = np.arccos(sigma)/3
hax4 = ax4.hist2d(theta, rho, bins=50, normed=True,\
       range=np.array([(0,np.pi/3),(0,2)]),\
        cmap='YlOrRd')
ax4.set_xlim([0,np.pi/3])
ax4.set_xticks([0,np.pi/6,np.pi/3])
ax4.set_xticklabels([r'$0$',r'$\pi/6$',r'$\pi/3$'])
ax4.set_ylim([0,2])
ax4.set_xlabel(r'$\theta$',labelpad=0)
ax4.set_ylabel(r'$\rho$',labelpad=0)
# ax4.legend([r'T='+str(T)])
axColorbar = inset_axes(ax4, width='3%',  height='30%', loc=1)
fig.colorbar(hax4[3], cax=axColorbar, ticks=[])

plt.show()
