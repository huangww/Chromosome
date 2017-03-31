import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import theory as th
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def gyrationRadius(r):
    """ compute the gyration radius of a configration

    :r: the position configuration
    :returns: float value of rg

    """
    rcm = r.mean(axis=0)
    rg2 = ((r - rcm)**2).sum()/len(r)
    return np.sqrt(rg2)

def posToRg():
    """ Transfer the data of position configuraton to the gyration radius
    :returns: null

    """
    fname = 'BeadRod_r_N100_F%g_T1_0.dat'% 0.1
    pos = np.loadtxt(fname).reshape([-1,100,3])
    rg = [gyrationRadius(r) for r in pos]
    fname = fname.replace('_r_', '_rg_')
    np.savetxt(fname, rg)
    return

def plotShape():
    """TODO: plot the shape profile and its marginal distribution
    :returns: null

    """
    fig = plt.figure(0,figsize=(8,4.2))
    F = 1
    fname = 'BeadRod_r_N100_F%g_T1_0.dat'%F
    data = np.loadtxt(fname)
    x = data[:,0]
    y = data[:,1]

    nullfmt = NullFormatter()         # no labels

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    ax = fig.add_axes(rect_scatter)
    axHistx = fig.add_axes(rect_histx)
    axHisty = fig.add_axes(rect_histy)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHistx.yaxis.set_major_formatter(nullfmt)
    axHisty.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    hax = ax.hist2d(x, y, bins=100, range=[[-10, 50], [-15, 15]], normed=True, cmap = 'Blues', vmin=0, vmax=0.02)
    axins = inset_axes(ax, width="30%",  height="5%", loc=2)
    fig.colorbar(hax[3], cax=axins, orientation='horizontal',ticks=[])

    # now determine nice limits by hand:
    ax.set_xlim((-10, 50))
    ax.set_ylim((-15, 15))
    ax.set_xlabel(r'$z$', fontsize=20, labelpad=0)
    ax.set_ylabel(r'$x$', fontsize=20, labelpad=0)


    xbins = np.linspace(x.min(), x.max(), 50)
    axHistx.hist(x, bins=xbins, normed=True)
    ybins = np.linspace(y.min(), y.max(), 50)
    axHisty.hist(y, bins=ybins, normed=True,orientation='horizontal')

    # now plot the theoretical results for the histogram
    axHistx.plot(xbins, th.particle_density_z(xbins, 100, 1/F),'r-',lw=2)
    axHisty.plot(th.particle_density_x(ybins, 100, 1/F), ybins, 'r-',lw=2)

    axHistx.set_xlim(ax.get_xlim())
    axHisty.set_ylim(ax.get_ylim())
    axHistx.set_ylabel(r'$p(z)$', fontsize=20, labelpad=0)
    axHisty.set_xlabel(r'$p(x)$', fontsize=20, labelpad=0)

    plt.show()
    return

def plotPdfRg():
    """ plot the distribution of gyration radius
    :returns: null

    """
    plt.rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 15 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    fig = plt.figure(0, figsize=(5, 3.6))
    fig.subplots_adjust(left=0.13, right =0.95,\
            bottom=0.13, top =0.95, wspace=0.25)
    ax = fig.add_subplot(111)
    force = [0.001, 0.03, 0.06, 0.1]
    for F in force:
        fname = 'BeadRod_rg_N100_F%g_T1_0.dat'%F
        rg = np.loadtxt(fname)
        labelString = r'$F=%g$'% F
        ax.hist(rg, bins=50, normed=True, alpha=0.25, label=labelString)
    ax.set_xlabel(r'$R_g$')
    ax.set_ylabel(r'$p(R_g)$',labelpad=0)
    ax.set_ylim([0,1.2])
    ax.legend(ncol=2, fontsize=15)
    plt.show()
    

if __name__ == "__main__":
    # plotShape()
    # posToRg()
    plotPdfRg()
