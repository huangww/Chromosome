import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx
from theory import *

fig = plt.figure(0,figsize=(10,4))
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 12 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text',usetex = True)

fig.subplots_adjust(left=0.1, right =0.95,\
        bottom=0.15, top =0.95, wspace=0.25)

sp1 = fig.add_subplot(121)
sp2 = fig.add_subplot(122)

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
cax = fig.add_axes([0.42, 0.65, 0.02, 0.25])
fig.text(0.425,0.6, r"$\tilde{T}$")
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
sp1.set_ylabel(r"$\left<z_i\right>/a$",fontsize=14)

# subplot 2
pos = [10,50,100,150]
Tlog = np.logspace(-2,3,1000)
cMap = plt.get_cmap('winter')
cNorm = colors.Normalize(vmin = 0, vmax = max(pos))
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)
cset = ['Cyan','SkyBlue','DodgerBlue','Blue']
vard = np.zeros(len(Tlog))
for s,c in zip(pos, cset):
    count = 0
    for T in Tlog:
        varr = r_var(index, N, 2*T)
        vard[count] = varr[s]
        count = count + 1
    colorVar = scalarMap.to_rgba(s)
    sp2.plot(Tlog,np.sqrt(2*vard), color=c)

sp2.fill_between(Tlog, 0, 2, \
        facecolor = 'Blue', alpha = 0.2)
# sp2.plot([0.5,0.5],[0,2],'k--')
fig.text(0.88,0.35, r"$i=10$")
fig.text(0.88,0.62, r"$i=50$")
fig.text(0.88,0.76, r"$i=100$")
fig.text(0.88,0.86, r"$i=150$")
sp2.set_xlim([0.01,1000])
# sp2.set_ylim([0,10])
sp2.set_xscale('log')
sp2.set_xlabel(r'$\tilde{T}$')
sp2.set_ylabel(r"$var[\mathbf{d}]^{1/2}/a$",fontsize=14)

fig.savefig('figure3.pdf')
plt.show()
