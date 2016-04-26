import numpy as np
import matplotlib.pyplot as plt

dataDir = 'data/'

def Plot1DTauScaling():
    NList = [50, 100, 200, 500]
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    for N in NList:
        fname = dataDir + 'tauTemp1D_N%g.dat'%N
        data = np.loadtxt(fname)
        # x = data[:,0]/(N**(1.5))
        x = data[:,0]
        y = data[:,1]/(N**2)
        ax.plot(x, y, '-o')
    ax.set_xlabel(r'$N^{-3/2}/F$', fontsize=16, labelpad=1)
    ax.set_ylabel(r'$<\tau>/N^2$', fontsize=16, labelpad=1)
    ax.ticklabel_format(fontsize=16)
    ax.legend(['50','100','200','500'])
    ax.set_xscale('log')
    # ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    # ax.set_xlim(0, 2)
    plt.show()
   
def PlotTauScaling():
    NList = [50, 100, 200, 500]
    fig = plt.figure(0)
    ax = fig.add_subplot(121)
    for N in NList:
        fname = dataDir + 'tauTemp1D_N%g.dat'%N
        data = np.loadtxt(fname)
        x = data[:,0]/(N**(1.5))
        y = data[:,1]/(N**2)
        ax.plot(x, y, '-o')
    x = np.linspace(0.01, 10.0, 100)
    ax.plot(x, np.exp(x/30.)-1)
    ax.set_xlabel(r'$N^{-3/2}/F$', fontsize=16, labelpad=1)
    ax.set_ylabel(r'$<\tau>/N^2$', fontsize=16, labelpad=1)
    ax.ticklabel_format(fontsize=16)
    # ax.legend(['50','100','200','500'])
    # ax.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
    ax.set_yscale('log')
    # ax.set_xlim(0, 10)
    ax = fig.add_subplot(122)
    data = np.loadtxt('data/tauForce3D.dat')
    ax.plot(data[:,0], data[:,1], '-o')
    ax.set_xlabel(r'$1/F$', fontsize=16, labelpad=1)
    ax.set_ylabel(r'$<\tau>$', fontsize=16, labelpad=1)
    ax.ticklabel_format(fontsize=16)
    plt.show()
   

def main():
    Plot1DTauScaling()
    # PlotTauScaling()

if __name__ == "__main__":
    main()
    
    
