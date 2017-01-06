import numpy as np
import matplotlib.pyplot as plt
import statsmodels.tsa.stattools as ss
from scipy.signal import savgol_filter
import itertools as it
import os 
import glob
import sys


dataDir = 'data/'

def ACF(x):
    acf= ss.acf(x, nlags=len(x), fft=True)
    end = next((i for i,v in enumerate(acf) if v<1e-6), \
            len(x)-1)
    return acf[:end]

def MergeACF(fname):
    # fname = dataDir + 'rg_N100_T'+str(T)+ \
    #         '_'+str(i)+'.dat'
    fileList = glob.glob(fname)
    tlist = []
    acflist = []
    for f in fileList:
        data = np.loadtxt(f)
        t = data[:,0]
        x = data[:,1]
        acf = ss.acf(x, nlags=len(x), fft=True)
        tlist.append(t)
        acflist.append(acf)
    arrlen = [len(x) for x in acflist]
    minlen = np.array(arrlen).min()
    acf = [x[:minlen] for x in acflist]
    acf = np.array(acf).mean(axis=0)
    end = next(i for i,v in enumerate(acf) if v<1e-6)
    acf = acf[:end]
    t = np.array(tlist[0][:end])
    return t, acf

def FindFitRange(acf):
    interval = len(acf)/300
    inter = interval if interval > 1 else 1
    x = np.log(acf[::inter])
    dx1 = np.maximum(0, savgol_filter(x, 5, 2, deriv=1))
    dx1sum = [dx1[:i+1].sum() for i in range(len(dx1))]
    end = next((i for i,v in enumerate(dx1sum) if v>0.1), \
            len(dx1sum)-1)
    acf = acf[:max(end*inter,30)]
    interval = len(acf)/300
    inter = interval if interval > 1 else 1
    x = np.log(acf[::inter])
    dx2abs = np.abs(savgol_filter(x, 5, 2, deriv=2))
    cutoff = len(x)/5
    fitLen = len(x)/3
    d2sum = [dx2abs[i:i+fitLen].sum() for i in \
            range(cutoff, len(x)-cutoff-fitLen)]
    t0 = (np.array(d2sum).argmin() + cutoff) * inter
    t1 = t0 + fitLen * inter
    return t0, t1

def FitTau(t, acf):
    t0, t1 = FindFitRange(acf)
    tFit = t[t0:t1]
    acfFit = acf[t0:t1]
    coefs = np.polyfit(tFit, np.log(acfFit), 1)
    tau = -1./coefs[0]
    if tau > 1e6 or tau < 0:
        tau = np.nan
    print 'tau =', tau
    return tau

def PlotFig(t, acf, save=False, figName='fig/acfFit.pdf'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(t, acf,'-')
    interval = len(acf)/300
    inter = interval if interval > 1 else 1
    x = np.log(acf[::inter])
    dx1 = np.maximum(0, savgol_filter(x, 5, 2, deriv=1))
    dx1sum = [dx1[:i+1].sum() for i in range(len(dx1))]
    ax.plot(t[::inter], dx1sum, 'r-', label='dx1_sum')
    t0, t1 = FindFitRange(acf)
    tFit = t[t0:t1]
    acfFit = acf[t0:t1]
    coefs = np.polyfit(tFit, np.log(acfFit), 1)
    tau = -1./coefs[0]
    ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]))
    ax.fill_betweenx(acf, t[t0], t[t1], alpha=0.1)
    ax.set_yscale('log')
    ax.set_xlabel(r'$\delta t$')
    ax.set_ylabel(r'$ACF$')
    ax.set_title('tau = %g'%tau)
    ax.set_ylim(1e-3, 1)
    ax.set_xlim(0, min(3*t[len(tFit)*3]/2, t[-1]))
    ax.axvline(t[len(tFit)*3])
    if save: fig.savefig(figName)
    plt.show()
    plt.close(fig)

def GetTau(t, x, saveFig=False, figName='fig/acfFit.pdf'):
    acf = ACF(x)
    t = t[:len(acf)]
    tau = FitTau(t, acf)
    if saveFig: PlotFig(t, acf, saveFig, figName)
    return tau

def ShowDataACF(data):
    t = data[:,0]
    x = data[:,1]
    acf = ACF(x)
    t = t[:len(acf)]
    PlotFig(t, acf)

def Demo():
    plt.close('all')
    fig = plt.figure(0)

    T = 1
    fname = 'data/acf1D_N100_T'+str(T)+'_c1.dat'
    t, acf = np.loadtxt(fname)
    ax = fig.add_subplot(131)
    ax.plot(t, acf,'-')
    t0, t1 = FindFitRange(acf)
    tFit = t[t0:t1]
    acfFit = acf[t0:t1]
    coefs = np.polyfit(tFit, np.log(acfFit), 1)
    ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]), 'r-')
    ax.set_yscale('log')
    ax.set_xlabel(r'$t$',fontsize=15, labelpad=1)
    ax.set_ylabel(r'$ACF$',fontsize=15, labelpad=1)
    ax.set_ylim(1e-3, 1)
    ax.set_title(r'$1/F=1$')

    T = 100
    fname = 'data/acf1D_N100_T'+str(T)+'_c1.dat'
    t, acf = np.loadtxt(fname)
    ax = fig.add_subplot(132)
    ax.plot(t, acf,'-')
    t0, t1 = FindFitRange(acf)
    tFit = t[t0:t1]
    acfFit = acf[t0:t1]
    coefs = np.polyfit(tFit, np.log(acfFit), 1)
    ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]), 'r-')
    ax.set_yscale('log')
    ax.set_xlabel(r'$t$', fontsize=15, labelpad=1)
    ax.set_ylabel(r'$ACF$', fontsize=15, labelpad=1)
    ax.set_ylim(1e-3, 1)
    ax.set_title(r'$1/F=100$')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

    T = 1000
    fname = 'data/acf1D_N100_T'+str(T)+'_c1.dat'
    t, acf = np.loadtxt(fname)
    ax = fig.add_subplot(133)
    ax.plot(t, acf,'-')
    t0, t1 = FindFitRange(acf)
    tFit = t[t0:t1]
    acfFit = acf[t0:t1]
    coefs = np.polyfit(tFit, np.log(acfFit), 1)
    ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]), 'r-')
    ax.set_yscale('log')
    ax.set_xlabel(r'$t$', fontsize=15, labelpad=1)
    ax.set_ylabel(r'$ACF$', fontsize=15, labelpad=1)
    ax.set_ylim(1e-3, 1)
    ax.set_title(r'$1/F=1000$')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.show()

def AddToFile(fname, T, tau):
    f = open(fname, 'ab')
    np.savetxt(f, (T, tau), newline=' ')
    f.write('\n')
    f.close()

def GetParaList():
    fname = dataDir + 'rg1D_N500_*.dat'
    fnameList = glob.glob(fname)
    NList = []
    TList = []
    for fname in fnameList:
        paraString = fname.replace('.','_').split('_')
        N = int(paraString[1][1:])
        T = float(paraString[2][1:])
        NList.append(N)
        TList.append(T)
    return np.unique(NList), np.unique(TList)


def main():
    # Tarr = list(np.linspace(0.2, 10, 99)) + \
    #          list(np.linspace(15, 100, 18))
    # Narr = [100]
    Tarr = [10]
    Narr = [100]
    # Narr, Tarr = GetParaList()
    for (N, T) in it.product(Narr, Tarr):
        # fname = dataDir+'rg1D_N%g_T%g_*.dat'%(N,T)
        # t, acf = MergeACF(fname)
        fname = dataDir+'rg1D_N%g_T%g.dat'%(N,T)
        data = np.loadtxt(fname)
        t, acf = MergeACF(data)
        # fname = 'data/acf1D_N100_T'+str(T)+'_c'+str(c)+'.dat'
        # np.savetxt(fname, (t, acf))
        # t, acf = np.loadtxt(fname)
        print N, T, len(acf)
        tau = FitTau(t, acf)
        PlotFig(t, acf)
        fname = 'data/tauForce_N'+str(N)+'.dat'
        # AddToFile(fname, T, tau)
    # Demo()
     
if __name__ == "__main__":
    main()
