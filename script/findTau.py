#!/usr/bin/env python
# encoding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from theory import *
from jsd import JSD
from scipy.interpolate import splrep, splev

def GetWatchPara(posMean, i0=1, mtd='posI0'):
    if mtd == 'posI0':
        watchPara = posMean[:,i0,0]
        zi0 = watchPara[-len(watchPara)/10:].mean()
        watchPara = np.fabs(watchPara - zi0)
        t = np.arange(len(watchPara)) 
    elif mtd == 'jsd':
        Q = z_mean(N, T)
        Q = np.maximum(Q,1e-6)
        t = np.arange(len(posMean))
        watchPara = np.zeros(len(posMean))
        for i in range(len(t)):
            P = posMean[i,:,0]
            P = np.maximum(P,1e-6)
            watchPara[i] = JSD(P, Q)
        watchPara = watchPara/watchPara.max()
    return t, watchPara

def FindLinear(x, y):
    f = splrep(x, y)
    d2 = splev(x, f, der=2)
    d2 = np.fabs(d2 / d2.max())
    end = next(i for i,v in enumerate(d2) if v>0.1)
    start = int(end/2)
    return x[start:end], y[start:end]

def FindTau(posMean):
    t, watchPara = GetWatchPara(posMean)
    x, y = FindLinear(t, watchPara)
    coefs = np.polyfit(x,y,1)
    tau = -1.0/coefs[0]
    return tau

def main():
    N = 100  # number of bead
    T = 1
    figDir = 'fig/'
    dataDir = 'data/'

    fname = dataDir+'mean_N'+str(N)+'_T'+str(T)+'.dat'
    posMean = np.loadtxt(fname).reshape(-1,N,3)

    fig = plt.figure(0, figsize=(9,6))
    ax = fig.add_subplot(111)

    t, watchPara = GetWatchPara(posMean, i0=50)
    markers, = ax.plot(t, watchPara, 'o')

    x, y = FindLinear(t, np.log(watchPara))
    print len(x)
    coefs = np.polyfit(x,y,1)
    yFit = t * coefs[0] + coefs[1]
    tau = -1.0/coefs[0]
    labelString =  r'$\tau=$' + str(tau)
    ax.plot(t, np.exp(yFit), color=markers.get_color(), label=labelString)

    # fname = dataDir + 'tau_N'+str(N)+'_T'+str(T)+'.dat'
    # np.savetxt(fname, np.vstack((t.T, watchPara.T)).T)

    ax.set_xlabel(r'$t$')
    ax.set_yscale('log')
    ax.set_ylabel(r"$d$")
    ax.legend()
    fname = figDir + 'tau_N'+str(N)+'_T'+str(T)+'.pdf'
    plt.savefig(fname)
    plt.close(fig)

if __name__ == '__main__':
    main()

