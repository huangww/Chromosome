#!/usr/bin/env python
# encoding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import glob
from theory import *
from jsd import JSD
from scipy.interpolate import splrep, splev
font = {'family' : 'scans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 12}
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text', usetex = True)

def GetWatchPara(data, i0=1, mtd='posI0'):
    t = data[:, 0]
    if mtd == 'posI0':
        watchPara = data[:,i0+1]
        zi0 = watchPara[-len(watchPara)/100:].mean()
        # zi0 = mean_zi_1D(i0, len(data[0])-1, 1)
        # zi0 = var_zi_1D(i0, len(data[0])-1, 1)
        watchPara = np.fabs(watchPara - zi0)
    elif mtd == 'jsd':
        posMean = data[:,1:]
        Q = z_mean_1D(N, T)
        Q = np.maximum(Q,1e-6)
        watchPara = np.zeros(len(t))
        for i in range(len(t)):
            P = posMean[i,:]
            P = np.maximum(P,1e-6)
            watchPara[i] = JSD(P, Q)
    return t, watchPara

def FindLinear(i0, x, y):
    f = splrep(x, y)
    d2 = splev(x, f, der=2)
    d2 = np.fabs(d2 / d2.max())
    end = next(i for i,v in enumerate(d2) if v>0.1)
    start = int(i0*9)
    # start = int(end/2)
    return x[start:end], y[start:end]

def FindTau(data):
    t, watchPara = GetWatchPara(data)
    x, y = FindLinear(10, t, np.log(watchPara))
    coefs = np.polyfit(x,y,1)
    tau = -1.0/coefs[0]
    return tau


if __name__ == '__main__':
    N = 100  # number of bead
    T = 10
    figDir = 'fig/'
    dataDir = 'data/'

    fname = dataDir + 'asep_mean_N' + str(N) \
            + '_T' + str(T) + '.dat'
    data = np.loadtxt(fname)

    fig = plt.figure(0, figsize=(9,6))
    ax = fig.add_subplot(111)

    t, watchPara = GetWatchPara(data, i0=50)
    markers, = ax.plot(t, watchPara, 'o')

    x, y = FindLinear(10, t, np.log(watchPara))
    print len(x)
    coefs = np.polyfit(x,y,1)
    yFit = t * coefs[0] + coefs[1]
    tau = -1.0/coefs[0]
    labelString =  r'$\tau=$' + str(tau)
    ax.plot(t, np.exp(yFit), color=markers.get_color(), label=labelString)

    # fname = 'tauN'+str(N)+'T'+str(T)+'I'+str(i0)+'.dat'
    # np.savetxt(figDir + fname, np.vstack((t.T, watchPara.T)).T)

    ax.set_xlabel(r'$t$')
    ax.set_yscale('log')
    ax.set_ylabel(r"$d$")
    ax.legend()
    fname = figDir + 'tau_1D_N'+str(N)+'_T'+str(T)+'.pdf'
    plt.savefig(fname)
    plt.show()
