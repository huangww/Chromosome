import numpy as np
import matplotlib.pyplot as plt
import statsmodels.tsa.stattools as ss
from scipy.interpolate import splrep, splev
import itertools as it
import os 
import glob


dataDir = 'data/'

def ACF(x):
    acf= ss.acf(x, nlags=len(x), fft=True)
    end = next((i for i,v in enumerate(acf) if v<0), len(x)-1)
    return acf[:end]

def GetACF0(data):
    acf = ACF(data[:,1])
    t = data[:len(acf),0]
    return t, acf

def GetACF(T, c):
    fname = dataDir + 'rg1D_N500_T'+str(T)+ \
            '_c'+str(c)+'_*.dat'
    print fname
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
    end = next(i for i,v in enumerate(acf) if v<0)
    acf = acf[:end]
    t = np.array(tlist[0][:end])
    return t, acf

def FindFitRange(acf):
    t0 = len(acf)/10
    t = np.arange(len(acf))
    f = splrep(t, np.log(acf))
    d1 = splev(t, f, der=1)
    d2 = splev(t, f, der=2)
    d2 = np.fabs(d2)
    t11 = next((i for i,v in enumerate(d1[t0:]) if v>3e-3),
            len(acf)-1)
    t12 = next((i for i,v in enumerate(d2[t0:]) if v>3e-3), 
            len(acf)-1)
    d2sum = [sum(d2[t0:i+1]) for i in range(t0,len(d2))]
    t13 = next((i for i,v in enumerate(d2sum) if v>1e-1), 
            len(acf)-t0-1)
    t1 = min(t11, t12, t13) + t0
    t0 = max(d2[:t1/3].argmin(), len(acf)/10)
    return t0, t1

def FitTau(t, acf):
    t0, t1 = FindFitRange(acf)
    tFit = t[t0:t1]
    acfFit = acf[t0:t1]
    coefs = np.polyfit(tFit, np.log(acfFit), 1)
    tau = -1./coefs[0]
    print 'tau =', tau
    return tau

def PlotFig(t, acf):
    plt.close('all')
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    ax.plot(t, acf,'-')
    t0, t1 = FindFitRange(acf)
    # ax.fill_betweenx(acf, t[t0], t[t1], alpha=0.1)
    tFit = t[t0:t1]
    acfFit = acf[t0:t1]
    coefs = np.polyfit(tFit, np.log(acfFit), 1)
    tau = -1./coefs[0]
    ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]))
    ax.set_yscale('log')
    ax.set_xlabel(r'$\delta t$')
    ax.set_ylabel(r'$ACF$')
    ax.set_ylim(1e-3, 1)
    ax.set_title(r'$1/F=1$')
    # f = splrep(np.arange(len(acf)), np.log(acf))
    # d1 = splev(np.arange(len(acf)), f, der=1)
    # d2 = splev(np.arange(len(acf)), f, der=2)
    # d2 = np.fabs(d2)
    # d2sum = [sum(d2[t0:i+1]) for i in range(t0,len(d2))]
    # ax.plot(t, d2)
    # ax.plot(t, d1)
    # ax.plot(t[t0:], d2sum)
    # ax.axhline(3e-3)
    # ax.axhline(1e-1)
    # ax.text(2*t[-1]/3, 0.5, r'$\tau=$'+str(int(tau)))
    plt.show()
    fig.savefig('fig/acfFit.pdf')

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
    ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]))
    ax.set_yscale('log')
    ax.set_xlabel(r'$\delta t$')
    ax.set_ylabel(r'$ACF$')
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
    ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]))
    ax.set_yscale('log')
    ax.set_xlabel(r'$\delta t$')
    ax.set_ylabel(r'$ACF$')
    ax.set_ylim(1e-3, 1)
    ax.set_title(r'$1/F=100$')

    T = 1000
    fname = 'data/acf1D_N100_T'+str(T)+'_c1.dat'
    t, acf = np.loadtxt(fname)
    ax = fig.add_subplot(133)
    ax.plot(t, acf,'-')
    t0, t1 = FindFitRange(acf)
    tFit = t[t0:t1]
    acfFit = acf[t0:t1]
    coefs = np.polyfit(tFit, np.log(acfFit), 1)
    ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]))
    ax.set_yscale('log')
    ax.set_xlabel(r'$\delta t$')
    ax.set_ylabel(r'$ACF$')
    ax.set_ylim(1e-3, 1)
    ax.set_title(r'$1/F=1000$')

    plt.show()

def AddToFile(fname, T, tau):
    f = open(fname, 'ab')
    np.savetxt(f, (T, tau), newline=' ')
    f.write('\n')
    f.close()

def main():
    Tarr = np.linspace(0.2, 10, 99) 
    Narr = [100, 500, 1000]
    for (N, T) in it.product(Narr, Tarr):
        fname = 'data/rg1D_N%g_T%g.dat'%(N,T)
        data = np.loadtxt(fname)
        # t, acf = GetACF(T, c)
        t, acf = GetACF0(data)
        # fname = 'data/acf1D_N100_T'+str(T)+'_c'+str(c)+'.dat'
        # np.savetxt(fname, (t, acf))
        # t, acf = np.loadtxt(fname)
        print N, T 
        tau = FitTau(t, acf)
        # PlotFig(t, acf)
        fname = 'data/tauStongForceN'+str(N)+'.dat'
        AddToFile(fname, T, tau)
    # Demo()
     
if __name__ == "__main__":
    main()
