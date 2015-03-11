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

sp = fig.add_subplot(111)

N = 300
index = np.arange(0,N)
dataDir = './'


Teff = [2,20]
# cMap = plt.get_cmap('gist_heat_r')
cMap = plt.get_cmap('autumn_r')
cNorm = colors.Normalize(vmin = 0, vmax = max(Teff))
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)
for T in Teff:
    fileName = dataDir + 'MD_N'+str(N)+'_T'+str(T)+'_xyz_varxyz.dat'
    data = np.loadtxt(fileName)
    colorVar = scalarMap.to_rgba(T)
    line, = sp.plot(index,z_var(index,N,int(T)),color=colorVar)
    linez, = sp.plot(index[::10],data[::10,3],'o',color=colorVar)
    line, = sp.plot(index,xy_var(index,N,int(T))/2.,'--',color=colorVar)
    linex, = sp.plot(index[::10],data[::10,4],'*',color=colorVar)
    liney, = sp.plot(index[::10],data[::10,5],'^',color=colorVar)
    # plot colorbar legend
cax = fig.add_axes([0.83, 0.66, 0.04, 0.25])
fig.text(0.835,0.6, r"$\tilde{T}$")
cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
cb.set_ticks([0,10,20])
cb.set_ticklabels([r'$0$',r'$5$',r'$10$'])
# cb.set_ticklabels([0,25,50])
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=2.0, headwidth=8.0))

sp.set_xlim([0,300])
sp.set_ylim([0,20])
sp.set_xlabel(r'Bead index $i$')
sp.set_ylabel(r"$var[x_i,y_i,z_i]/a^2$")
sp.legend((linex, liney, linez), (r'$x$',r'$y$',r'$z$'), loc = 'upper left', handlelength=1)

fig.savefig('supfig.pdf')
plt.show()
