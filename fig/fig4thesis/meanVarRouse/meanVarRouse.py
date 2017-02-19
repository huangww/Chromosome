import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.colorbar as colorbar
import matplotlib.cm as cmx
import theory as th


N = 100

def plotMeanVar():
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
    Teff = [5,10,20]
    # cMap = plt.get_cmap('gist_heat_r')
    cMap = plt.get_cmap('winter')
    cNorm = colors.Normalize(vmin = 0, vmax = max(Teff))
    scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cMap)
    
    for T in Teff:
        colorVar = scalarMap.to_rgba(T)
        fname = 'BeadSpring_meanVar_N100_T%g.dat'%T
        data = np.loadtxt(fname)
        zmean = data[0]
        rvar = data[3] + data[4] + data[5]
        zmean, rvar = (np.append(zmean, zmean[0]), np.append(rvar, rvar[0]))
        ax1.plot(index[::5], zmean[::5], 'o', color = colorVar)
        ax2.plot(index[::5], rvar[::5], 'o', color = colorVar)
        rmean = th.rouse_mean(N, T)
        rmean = np.append(rmean, rmean[0])
        lineSpring, = ax1.plot(index, rmean, color=colorVar)
        zmeanRod = [th.mean_zi(i, N, T) for i in index]
        lineRod, = ax1.plot(index, zmeanRod, '--', color = colorVar) 
        ax2.plot(index, index*(N-index)/float(N), color = colorVar)

    # plot line of Tinf
    ax1.plot([0,N],[0,0],'k--')
    # zvar = [var_zi_1D(i, N, 1000000) for i in index]
    rvar = index*(N-index)/float(N)
    ax2.plot(index, rvar,'k--')

    # add colorbar legend
    cax = fig.add_axes([0.42, 0.66, 0.02, 0.25])
    fig.text(0.42,0.6, r"$\tilde{T}$", fontsize=15)
    # fig.text(0.405,0.59, r"$\frac{k_B T}{2Fa}$", fontsize=15)
    cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
    cb.set_ticks([0,10,20])
    for T in Teff:
        colorVar = scalarMap.to_rgba(T)
        cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=0., headwidth=6.0))

    cax = fig.add_axes([0.892, 0.66, 0.02, 0.25])
    fig.text(0.892,0.6, r"$\tilde{T}$",fontsize=15)
    # fig.text(0.88,0.59, r"$\frac{k_B T}{2Fa}$", fontsize=15)
    cb = colorbar.ColorbarBase(cax, cmap = cMap, norm = cNorm)
    cb.set_ticks([0,10,20])
    for T in Teff:
        colorVar = scalarMap.to_rgba(T)
        cax.annotate('', xy=(-0.0, T/float(max(Teff))), xytext=(-1.0, T/float(max(Teff))), arrowprops=dict(facecolor=colorVar,edgecolor='none',width=0.0, headwidth=6.0))

    # set labels and ticks
    ax1.set_xlabel(r'$\mathrm{Bead}\ \mathrm{index}\ i$')
    # ax1.set_xticks([0, 20, 40, 60])
    ax1.set_ylabel(r"$\left<z_i\right>/a$")
    # ax1.set_yticks([0, 10, 20, 30, 40])
    # ax1.set_ylim([-2, 30])
    ax2.set_xlabel(r'$\mathrm{Bead}\ \mathrm{index}\ i$')
    # ax2.set_xticks([0, 100, 200, 300])
    ax2.set_ylabel(r"$\mathrm{var}\left[\mathbf{r}_i\right]/a^2$")
    # ax2.set_yticks([0, 20, 40, 60, 80])
    # ax2.set_ylim([0, 10])
    # ax2.set_ylim([0, 50])

    # fig.set_tight_layout(True)
    fig.savefig('meanVarRouse.pdf')
    plt.show()

def main():
    plotMeanVar()

if __name__ == "__main__":
    main()
