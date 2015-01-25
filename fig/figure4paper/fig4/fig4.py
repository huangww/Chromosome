import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx
from theory import *

fig = plt.figure(0,figsize=(4,3))
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 12 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text',usetex = True)
fig.subplots_adjust(left=0.15, right =0.95,\
        bottom=0.15, top =0.95, wspace=0.25)
ax = fig.add_subplot(111)

N = 50
cm = 19     # Centromere location
dataDir = './'
index = np.linspace(1,N,N)
Teff = [1,10,50]

cMap = plt.get_cmap('autumn_r')
cNorm = colors.Normalize(vmin = 0, vmax = max(Teff))
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    varz = var_z(N, T)
    ax.plot(index,varz,'--', color=colorVar)
    varzcm = var_z_cm(cm, N, T)
    ax.plot(index, varzcm, color=colorVar)
    fileName = dataDir + 'MD_N98_T'+str(T) + \
            '_pseudo_xyz_varxyz.dat'
    data = np.loadtxt(fileName)
    ax.plot(index[::2],data[3,::2],'o', color=colorVar)

# plot High T limit
T = 5000
varz = var_z(N, T)
ax.plot(index,varz,'k--')

# add colorbar legend
cax = fig.add_axes([0.85, 0.66, 0.04, 0.25])
fig.text(0.855,0.6, r"$\tilde{T}$")
cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
cb.set_ticks([0,25,50])
cb.set_ticklabels([r'$0$',r'$25$',r'$50$'])
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=1.5, headwidth=6.0))

# ax.set_xlim([0,300])
# ax.set_ylim([0,120])
ax.set_xlabel(r'Bead index $i$')
ax.set_ylabel(r"$\left<z_i\right>/a$",fontsize=14)


fig.savefig('figure4.pdf')
plt.show()
