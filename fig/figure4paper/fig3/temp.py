import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx
from theory import *

fig = plt.figure(0,figsize=(5,4))
plt.rc('text',usetex = True)
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 16 }
plt.rc('lines', lw=2)
plt.rc('font', **font)

fig.subplots_adjust(left=0.15, right =0.95,\
        bottom=0.15, top =0.95, wspace=0.25)

sp1 = fig.add_subplot(111)

N = 300
index = np.arange(0,N)
dataDir = './'


# subplot 1
Teff = [10,20,50,100]
# cMap = plt.get_cmap('gist_heat_r')
cMap = plt.get_cmap('autumn_r')
cNorm = colors.Normalize(vmin = 0, vmax = max(Teff))
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)
for T in Teff:
    fileName = dataDir + 'MD_N'+str(N)+'_T'+str(T)+'_xyz_varxyz.dat'
    data = np.loadtxt(fileName)
    colorVar = scalarMap.to_rgba(T)
    line, = sp1.plot(index,z_mean(index,N,int(T)),color=colorVar)
    sp1.plot(index[::10],data[::10,0],'o',color=colorVar)
    sp1.fill_between(index, \
            z_mean(index,N,T)-np.sqrt(z_var(index,N,T)), \
            z_mean(index,N,T)+np.sqrt(z_var(index,N,T)), \
            facecolor = colorVar,\
            alpha = 0.2)
# plot colorbar legend
cax = fig.add_axes([0.82, 0.65, 0.04, 0.25])
fig.text(0.821,0.59, r"$\tilde{T}$")
cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
cb.set_ticks([0,50,100])
cb.set_ticklabels([r'$0$',r'$25$',r'$50$'])
# cb.set_ticklabels([0,25,50])
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=1.5, headwidth=6.0))

sp1.set_xlim([0,300])
sp1.set_ylim([0,120])
sp1.set_xlabel(r'Bead index $i$')
sp1.set_ylabel(r"$\left<z_i\right>/a$")

fig.savefig('figure3.pdf')
plt.show()
