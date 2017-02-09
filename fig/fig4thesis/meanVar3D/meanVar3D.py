import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx
from theory import *


N = 100

def meanVarData(fname):
    data = np.loadtxt(fname)
    pos = data.reshape([-1, N, 3])
    meanPos = pos.mean(axis=0)
    varPos = pos.var(axis=0)
    return meanPos, varPos

def plotMeanVarZ():
    fig = plt.figure(0, figsize=(10,4))
    plt.rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 15 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    fig.subplots_adjust(left=0.1, right =0.95,\
        bottom=0.18, top =0.95, wspace=0.25)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    
    index = np.arange(N+1)
    Teff = [10,20,50,100]
    # cMap = plt.get_cmap('gist_heat_r')
    cMap = plt.get_cmap('winter')
    cNorm = colors.Normalize(vmin = 0, vmax = max(Teff))
    scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)
    
    for T in Teff:
        fname = 'BeadRod_meanVar_N100_T%g.dat'%T
        data = np.loadtxt(fname)
        meanz, varz = (data[:,0], data[:,3])
        meanz, varz = (np.append(meanz, meanz[0]), np.append(varz, varz[0]))
        colorVar = scalarMap.to_rgba(T)
        zmean = [mean_zi(i, N, T) for i in index]
        zvar = [var_zi(i, N, T) for i in index]
        line, = ax1.plot(index, zmean, color = colorVar)
        ax1.plot(index[::5], meanz[::5], 'o', color = colorVar)
        ax2.plot(index, zvar, color = colorVar)
        ax2.plot(index[::5], varz[::5], 'o', color = colorVar)

    # plot line of Tinf
    ax1.plot([0,N],[0,0],'k--')
    # zvar = [var_zi_1D(i, N, 1000000) for i in index]
    zvar = [var_zi(i, N, 1000000) for i in index]
    ax2.plot(index, zvar,'k--')

    # add colorbar legend
    cax = fig.add_axes([0.42, 0.66, 0.02, 0.25])
    fig.text(0.42,0.6, r"$\tilde{T}$", fontsize=15)
    # fig.text(0.405,0.59, r"$\frac{k_B T}{2Fa}$", fontsize=15)
    cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
    cb.set_ticks([0,50,100])
    for T in Teff:
        colorVar = scalarMap.to_rgba(T)
        cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=0., headwidth=6.0))

    cax = fig.add_axes([0.892, 0.66, 0.02, 0.25])
    fig.text(0.892,0.6, r"$\tilde{T}$",fontsize=15)
    # fig.text(0.88,0.59, r"$\frac{k_B T}{2Fa}$", fontsize=15)
    cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
    cb.set_ticks([0,50,100])
    for T in Teff:
        colorVar = scalarMap.to_rgba(T)
        cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=0.0, headwidth=6.0))

    # set labels and ticks
    ax1.set_xlabel(r'$\mathrm{Bead}\ \mathrm{index}\ i$')
    # ax1.set_xticks([0, 20, 40, 60])
    ax1.set_ylabel(r"$\left<z_i\right>/a$")
    ax1.set_yticks([0, 10, 20, 30, 40])
    ax1.set_ylim([-2, 30])
    ax2.set_xlabel(r'$\mathrm{Bead}\ \mathrm{index}\ i$')
    # ax2.set_xticks([0, 100, 200, 300])
    ax2.set_ylabel(r"$\mathrm{var}\left[z_i\right]/a^2$")
    # ax2.set_yticks([0, 20, 40, 60, 80])
    ax2.set_ylim([0, 10])
    # ax2.set_ylim([0, 50])

    # fig.set_tight_layout(True)
    fig.savefig('meanVarZ3D.pdf')
    plt.show()

def plotVarXY():
    fig = plt.figure(0,figsize=(5,4))
    plt.rc('text',usetex = True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 15 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    fig.subplots_adjust(left=0.15, right =0.95,\
            bottom=0.18, top =0.95, wspace=0.25)
    ax = fig.add_subplot(111)

    index = np.arange(N+1)
    Teff = [1, 10, 20]
    # cMap = plt.get_cmap('gist_heat_r')
    cMap = plt.get_cmap('winter')
    cNorm = colors.Normalize(vmin = 0, vmax = max(Teff))
    scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)
    for T in Teff:
        fname = 'BeadRod_meanVar_N100_T%g.dat'%T
        data = np.loadtxt(fname)
        varx, vary, varz = (data[:,5], data[:,4], data[:,3])
        varx, vary, varz = (np.append(varx, varx[0]), np.append(vary, vary[0]), np.append(varz, varz[0]))
        zvar = [var_zi(i, N, T) for i in index]
        xyvar = [var_xyi(i, N, T)/2. for i in index]
        colorVar = scalarMap.to_rgba(T)
        line, = ax.plot(index, zvar, color=colorVar)
        linez, = ax.plot(index[::5],varz[::5],'o',color=colorVar)
        line, = ax.plot(index, xyvar,'--', color=colorVar)
        linex, = ax.plot(index[::5],varx[::5],'*',color=colorVar)
        liney, = ax.plot(index[::5],vary[::5],'^',color=colorVar)
    # plot colorbar legend
    cax = fig.add_axes([0.83, 0.66, 0.04, 0.25])
    fig.text(0.835,0.59, r"$\tilde{T}$")
    cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
    cb.set_ticks([0,10,20])
    # cb.set_ticklabels([r'$0$',r'$5$',r'$10$'])
    # cb.set_ticklabels([0,25,50])
    for T in Teff:
        colorVar = scalarMap.to_rgba(T)
        cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=0.0, headwidth=6.0))

    # ax.set_xticks([0,100,200,300])
    ax.set_ylim([0,10])
    ax.set_xlabel(r'$\mathrm{Bead}\ \mathrm{index}\ i$')
    ax.set_ylabel(r"$\mathrm{var}\left[x_i,y_i,z_i\right]/a^2$")
    ax.legend((linex, liney, linez), (r'$x$',r'$y$',r'$z$'), loc = 'upper left', handlelength=1, fontsize = 18)

    fig.savefig('xyVar3D.pdf')
    plt.show()
   

def main():
    # plotMeanVarZ()
    plotVarXY()

if __name__ == "__main__":
    main()
