import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx
from theory import *

fig = plt.figure(0,figsize=(5,4))
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 20 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text',usetex = True)
fig.subplots_adjust(left=0.15, right =0.95,\
        bottom=0.18, top =0.95, wspace=0.25)
ax = fig.add_subplot(111)

N = 300
cm = 100     # Centromere location
dataDir = './'
index = np.arange(N+1)
# Teff = [10,20,50,100]
Teff = [10,20,50]

cMap = plt.get_cmap('winter')
cNorm = colors.Normalize(vmin = 0, vmax = max(Teff))
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    varz = 2*var_z(N, T) # variance of distance
    varzcm = var_z_cm(cm, N, T)
    varz, varzcm = (np.append(varz, varz[0]), np.append(varzcm, varzcm[0]))
    # varz = np.concatenate((varz,[0]))
    ax.plot(index,varz,'--', color=colorVar)
    # varzcm = np.concatenate((varzcm,[0]))
    ax.plot(index, varzcm, color=colorVar)
    # fileName = dataDir + 'MD_N98_T'+str(T) + \
            # '_pseudo_xyz_varxyz.dat'
    fileName = dataDir + 'MD_N598_T'+str(T) + \
            '_cut_cm_xyz_varxyz.dat'
    data = np.loadtxt(fileName)
    ax.plot(index[::10],data[::10,3],'o', color=colorVar)

# plot High T limit
T = 10000
varz = 2*var_z(N, T)
varz = np.append(varz, varz[0])
ax.plot(index,varz,'k--')

# add colorbar legend
cax = fig.add_axes([0.84, 0.65, 0.04, 0.25])
fig.text(0.845,0.58, r"$\tilde{T}$")
cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
cb.set_ticks([0,20,50])
# cb.set_ticklabels([r'$0$',r'$10$',r'$25$'])
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=0.0, headwidth=6.0))

# ax.set_xticks([0,100, 200, 300])
# ax.set_ylim([0,120])
ax.set_xlabel(r'$\mathrm{Bead}\ \mathrm{index}\ i$')
ax.set_ylabel(r"$\mathrm{var}\left[d_{i,z}\right]/a^2$")

fig.savefig('figure4.pdf')
plt.show()
