import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx

def mean_n_1D(i, N, T):
    Teff = float(T)
    mu = (N-1)/2.0
    mean_n_1D = 1.0/(1+np.exp((i-mu)/Teff))
    return mean_n_1D

def var_n_1D(i, N, T):
    p = mean_n_1D(i, N, T)
    var_n_1D = p * (1 - p)
    return var_n_1D

def mean_zi_1D(i, N, T):
    Teff = float(T)
    mean_zi_1D = 0
    for i0 in range(i):
        mean_zi_1D += (2*mean_n_1D(i0, N, Teff)-1)
    return mean_zi_1D

def z_var_1D(m, n, N, T):
    z_var_1D = 0
    for i in range(m,n):
        z_var_1D += 4*var_n_1D(i, N, T)
    return z_var_1D

def var_zi_1D_add(i, N, T):
    Teff = float(T)
    var_zi_1D_add = 0
    for i0 in range(min(i, N-i)):
        var_zi_1D_add += 4*var_n_1D(i0, N, Teff)
    return var_zi_1D_add

def var_zi_1D(i, N, T):
    Teff = float(T)
    var_zi_1D = 0
    for i0 in range(i):
        v0s = z_var_1D(0,i0,N,T)
        vst = z_var_1D(i0,N,N,T)
        v0t = v0s + vst
        var_zi_1D = v0s*vst/v0t
    return var_zi_1D


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

N = 100
dataDir = './'
Teff = [10,20,30,50]
# index = np.arange(N)
index = np.arange(N+1)

# cMap = plt.get_cmap('gist_heat_r')
cMap = plt.get_cmap('winter')
cNorm = colors.Normalize(vmin = 0, vmax = max(Teff))
scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)

for T in Teff:
    fileName = dataDir + 'Asep_posMeanVar_N' + str(N) + \
            '_T' + str(T) + '.dat'
    meanz, varz = np.loadtxt(fileName)
    meanz, varz = (np.append(meanz, meanz[0]), np.append(varz, varz[0]))
    colorVar = scalarMap.to_rgba(T)
    zmean = [mean_zi_1D(i, N, T) for i in index]
    # zvar = [var_zi_1D(i, N, T) for i in index]
    zvar = [var_zi_1D(i, N, T) for i in index]
    line, = sp1.plot(index,zmean, color = colorVar)
    sp1.plot(index[::5], meanz[::5], 'o', color = colorVar)
    sp2.plot(index,zvar, color = colorVar)
    sp2.plot(index[::5], varz[::5], 'o', color = colorVar)

# plot line of Tinf
sp1.plot([0,N],[0,0],'k--')
# zvar = [var_zi_1D(i, N, 1000000) for i in index]
zvar = [var_zi_1D(i, N, 1000000) for i in index]
sp2.plot(index,zvar,'k--')

# add colorbar legend
cax = fig.add_axes([0.41, 0.66, 0.02, 0.25])
fig.text(0.41,0.59, r"$\tilde{T}$", fontsize=15)
# fig.text(0.405,0.59, r"$\frac{k_B T}{2Fa}$", fontsize=15)
cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
cb.set_ticks([0,25,50])
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=0., headwidth=6.0))

cax = fig.add_axes([0.885, 0.66, 0.02, 0.25])
fig.text(0.885,0.59, r"$\tilde{T}$",fontsize=15)
# fig.text(0.88,0.59, r"$\frac{k_B T}{2Fa}$", fontsize=15)
cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
cb.set_ticks([0,25,50])
for T in Teff:
    colorVar = scalarMap.to_rgba(T)
    cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=0.0, headwidth=6.0))

# set labels and ticks
sp1.set_xlabel(r'$\mathrm{Bead}\ \mathrm{index}\ i$')
# sp1.set_xticks([0, 20, 40, 60])
sp1.set_ylabel(r"$\left<z_i\right>/a$")
sp1.set_yticks([0, 10, 20, 30, 40])
sp1.set_ylim([-2, 40])
sp2.set_xlabel(r'$\mathrm{Bead}\ \mathrm{index}\ i$')
# sp2.set_xticks([0, 100, 200, 300])
sp2.set_ylabel(r"$\mathrm{var}\left[z_i\right]/a^2$")
# sp2.set_yticks([0, 20, 40, 60, 80])
sp2.set_ylim([0, 30])
# sp2.set_ylim([0, 50])

zvar = [var_zi_1D(i, N, 1000000) for i in index]
# fig.set_tight_layout(True)
fig.savefig('meanVar1D.pdf')
plt.show()
