import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx
from theory import *
   
fig = plt.figure(0,figsize=(10,4))
plt.rc('text', usetex=True)
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 20 }
plt.rc('lines', lw=2)
plt.rc('font', **font)

fig.subplots_adjust(left=0.1, right =0.95,\
        bottom=0.18, top =0.95, wspace=0.25)


sp1 = fig.add_subplot(121)
sp2 = fig.add_subplot(122)

N = 50
dataDir = './'
Teff = [1,5,10,50]
index = np.linspace(1,N,N)

# cMap = plt.get_cmap('gist_heat_r')
cMap = plt.get_cmap('autumn_r')
cNorm = colors.Normalize(vmin = 0, vmax = max(Teff))
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)

for T in Teff:
    fileName = dataDir + 'MD_N' + str(N) + \
        '_T' + str(T) + '_xyz_varxyz_pseudo_LJ.dat'
    data = np.loadtxt(fileName)
    posMean = data[:3]
    posVar = data[3:]
    colorVar = scalarMap.to_rgba(T)
    line, = sp1.plot(mean_z(N,T), color = colorVar)
    sp1.plot(posMean[0], 'o', color = colorVar)
    sp2.plot(var_z(N,T), color = colorVar)
    sp2.plot(posVar[0], 'o', color = colorVar)

# plot line of Tinf
sp1.plot([0,N],[0,0],'k--')
sp2.plot(var_z(N,10000),'k--')

# add colorbar legend
cax = fig.add_axes([0.41, 0.66, 0.02, 0.25])
fig.text(0.41,0.59, r"$\tilde{T}$")
cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
cb.set_ticks([0,20,50])
cb.set_ticklabels([r'$0$',r'$10$',r'$25$'])
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=1.5, headwidth=6.0))

cax = fig.add_axes([0.885, 0.66, 0.02, 0.25])
fig.text(0.885,0.59, r"$\tilde{T}$")
cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
cb.set_ticks([0,20,50])
cb.set_ticklabels([r'$0$',r'$10$',r'$25$'])
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=1.5, headwidth=6.0))

# set labels and ticks
sp1.set_xlabel(r'$\mathrm{Bead}\ \mathrm{index}\ i$')
# sp1.set_xticks([0, 100, 200, 300])
sp1.set_ylabel(r"$\left<z_i\right>/a$")
sp1.set_yticks([0, 5, 10, 15, 20, 25])
sp1.set_ylim([-1, 25])
sp2.set_xlabel(r'$\mathrm{Bead}\ \mathrm{index}\ i$')
# sp2.set_xticks([0, 100, 200, 300])
sp2.set_ylabel(r"$\mathrm{var}\left[z_i\right]/a^2$")
sp2.set_yticks([0, 2, 4, 6, 8, 10])

# fig.set_tight_layout(True)
fig.savefig('figLJ.pdf')
plt.show()
