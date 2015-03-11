import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx

def z_mean(i, N, T):
    z_mean = 2*T*np.log((1.0+np.exp(N/float(2*T)))/ \
            (np.exp(i/float(2*T))+np.exp((N-i)/float(2*T))))
    return z_mean

def z_var(i, N, T):
    z_var = 2*T*np.sinh((N-i)/float(2*T))*np.sinh(i/float(2*T))/ \
            (np.sinh(N/float(2*T))*np.cosh((N-2*i)/float(4*T))**2)
    return z_var

# def raw_var(i, N, T):
#     T = 2*T
#     raw_var = 2.0*T/(1.0+np.exp(2.0*(i-N/2)/float(T)))
#     return raw_var
#
# def z_var(i, N, T):
#     z_var = (raw_var(0,N,T)-raw_var(i,N,T))*(raw_var(i,N,T)-raw_var(N,N,T))/(raw_var(0,N,T)-raw_var(N,N,T))
#     return z_var
    
fig = plt.figure(0,figsize=(10,4))
plt.rc('text', usetex=True)
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 16 }
plt.rc('lines', lw=2)
plt.rc('font', **font)

fig.subplots_adjust(left=0.1, right =0.95,\
        bottom=0.15, top =0.95, wspace=0.25)


sp1 = fig.add_subplot(121)
sp2 = fig.add_subplot(122)

N = 300
dataDir = './'
Teff = [10,25,50,100]
index = np.linspace(1,N,N)

# cMap = plt.get_cmap('gist_heat_r')
cMap = plt.get_cmap('autumn_r')
cNorm = colors.Normalize(vmin = 0, vmax = max(Teff))
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)

for T in Teff:
    fileName = dataDir + 'visualization_T' + str(T) + \
            '_N' + str(N) + '.dat'
    data = np.loadtxt(fileName)
    meanz = data.mean(axis=0)
    varz = data.var(axis=0)
    colorVar = scalarMap.to_rgba(T)
    line, = sp1.plot(index,z_mean(index,N,T), color = colorVar)
    sp1.plot(index[::10], meanz, 'o', color = colorVar)
    sp2.plot(index,z_var(index,N,T), color = colorVar)
    sp2.plot(index[::10], varz, 'o', color = colorVar)

# plot line of Tinf
sp1.plot([0,N],[0,0],'k--')
sp2.plot(index,z_var(index,N,10000),'k--')

# add colorbar legend
cax = fig.add_axes([0.417, 0.65, 0.02, 0.25])
fig.text(0.42,0.59, r"$\tilde{T}$")
cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
cb.set_ticks([0,50,100])
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=1.5, headwidth=6.0))

cax = fig.add_axes([0.89, 0.65, 0.02, 0.25])
fig.text(0.893,0.59, r"$\tilde{T}$")
cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
cb.set_ticks([0,50,100])
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=1.5, headwidth=6.0))

# set labels and ticks
sp1.set_xlabel(r'Bead index $i$')
sp1.set_ylabel(r"$\left<z_i\right>/a$")
sp1.set_yticks([0, 40, 80, 120, 160])
sp1.set_ylim([-5, 160])
sp2.set_xlabel(r'Bead index $i$')
sp2.set_ylabel(r"$var[z_i]/a^2$")
sp2.set_yticks([0, 20, 40, 60, 80])

fig.savefig('figure2.pdf')
plt.show()
