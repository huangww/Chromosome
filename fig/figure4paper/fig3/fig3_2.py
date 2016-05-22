import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx
from theory import *

fig = plt.figure(0)
plt.rc('text',usetex = True)
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 20 }
plt.rc('lines', lw=2)
plt.rc('font', **font)

fig.subplots_adjust(left=0.1, right =0.95,\
        bottom=0.18, top =0.95, wspace=0.25)
ax = fig.add_subplot(111)

N = 300
index = np.arange(0,N)
dataDir = './'

Teff = [1,2,10,20,50,100,200]
pos = [10,50,100,150]
Tlog = np.logspace(-2,3,1000)
cMap = plt.get_cmap('winter')
cNorm = colors.Normalize(vmin = 0, vmax = max(pos))
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)
cset = ['LightSkyBlue','DodgerBlue','RoyalBlue','DarkBlue']
vard = np.zeros(len(Tlog))
ax.fill_between(Tlog, 0, 2, \
        facecolor = 'Blue', alpha = 0.1)
for s,c in zip(pos, cset):
    count = 0
    for T in Tlog:
        varr = r_var(index, N, 2*T)
        vard[count] = varr[s]
        count = count + 1
    colorVar = scalarMap.to_rgba(s)
    ax.plot(Tlog,np.sqrt(2*vard), color=c)
    ax.axhline(y=max(np.sqrt(2*vard)), linestyle='--', color=c)
    for T in Teff:
        fileName = dataDir + 'MD_N'+str(N)+'_T'+str(T)+'_xyz_varxyz.dat'
        data = np.loadtxt(fileName)
        ax.scatter(T/2., np.sqrt(2*sum(data[s,3:])), color=c)
    fileName = dataDir + 'MD_N'+str(N)+'_T02_xyz_varxyz.dat'
    data = np.loadtxt(fileName)
    ax.scatter(0.2/2., np.sqrt(2*sum(data[s,3:])), color=c)
    fileName = dataDir + 'MD_N'+str(N)+'_T05_xyz_varxyz.dat'
    data = np.loadtxt(fileName)
    ax.scatter(0.5/2., np.sqrt(2*sum(data[s,3:])), color=c)
    

# ax.plot([0.5,0.5],[0,2],'k--')
fig.text(0.8,0.34, r"$i=10$")
fig.text(0.8,0.575, r"$i=50$")
fig.text(0.8,0.70, r"$i=100$")
fig.text(0.8,0.82, r"$i=150$")
ax.set_xlim([0.01,1000])
ax.set_ylim([0,15])
ax.set_yticks([0,5,10,15])
ax.set_xscale('log')
ax.set_xlabel(r'$1/F$')
ax.set_ylabel(r"$\mathrm{var}\left[\mathbf{d}\right]^{1/2}/a$")

fig.savefig('figure3_1.pdf')
plt.show()
