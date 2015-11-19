import numpy as np
import matplotlib.pyplot as plt

dataDir = 'data/'
figDir = 'fig/'

def LoadData(N, T):
    fname = dataDir + 'rg1D_N' + str(N) + '_T' + '%g'%T + '.dat'
    data = np.loadtxt(fname)
    return data

def GetRg(data, norm=False):
    nbins = 500
    binrange = (0,1e5)
    binCount, t = np.histogram(data[:,0], bins=nbins, 
            range=binrange)
    rg, t = np.histogram(data[:,0], bins=nbins, 
            range=binrange, weights=data[:,1])
    rg = rg / binCount
    if norm:
        rgeq = rg[-len(rg)/10].mean()
        rg = np.fabs(rg - rgeq)
        rg = rg / rg.max()
        return t[:-1], rg
    else:
        return t[:-1], rg


def main():
    N = 100
    T = 10
    data = LoadData(N, T)
    t, rg = GetRg(data, norm=True)
    # t, rg = GetRg(data)
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    ax.plot(t, rg)

    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'$<R_g(t)>-<R_g^{eq}>$')
    plt.show()
    fname = figDir + 'rg1D_N' + str(N) + '_T' + '%g' %T + '.pdf'
    fig.savefig(fname)

if __name__ == "__main__":
    main()
