import numpy as np
import matplotlib.pyplot as plt
import statsmodels.tsa.stattools as ss
from scipy.interpolate import splrep, splev

Temp = [1, 30, 50, 100, 200, 1000]

def ACF(data):
    t = data[:, 0]
    x = data[:, 1]
    acf = ss.acf(x, nlags=len(x), fft=True)
    end = next(i for i,v in enumerate(acf) if v<0)
    return t[:end], acf[:end]

def ACFData():
    tList = []
    acfList = []
    for T in Temp:
        fname = 'data/rg1D_N100_T' + str(T) + '_eq.dat'
        data = np.loadtxt(fname)
        t, acf = ACF(data)
        tList.append(t)
        acfList.append(acf)
    return tList, acfList

def LoadAcfData():
    tList = []
    acfList = []
    for T in Temp:
        fname = 'data/acf_N100_T'+str(T)+'.dat'
        data = np.loadtxt(fname)
        t, acf = data[0], data[1]
        end = next((i for i,v in enumerate(acf) if v<0), 
                len(acf)-1)
        t, acf = t[:end], acf[:end]
        tList.append(t)
        acfList.append(acf)
    return tList, acfList

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
    t13 = next((i for i,v in enumerate(d2sum) if v>2e-1), 
            len(acf)-t0-1)
    t1 = min(t11, t12, t13) + t0
    t0 = max(d2[:t1/3].argmin(), len(acf)/10)
    return t0, t1

def PlotFig(t, acf, indicator=False):
    plt.close('all')
    fig = plt.figure(0)
    fig.subplots_adjust(left=0.05, right =0.95,\
        bottom=0.05, top =0.95, wspace=0.2)
    for i,T in enumerate(Temp):
        ax = fig.add_subplot(2, 3, i+1)
        t0, t1 = FindFitRange(acf[i])
        if T == 1:
            t1 = 300
        tFit = t[i][t0:t1]
        acfFit = acf[i][t0:t1]
        coefs = np.polyfit(tFit, np.log(acfFit), 1)
        tau = -1./coefs[0]
        ax.plot(t[i], acf[i])
        ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]),lw=2)
        if indicator:
            x = np.arange(len(acf[i]))
            f = splrep(x, np.log(acf[i]))
            d1 = splev(x, f, der=1)
            d2 = splev(x, f, der=2)
            d2 = np.fabs(d2)
            d2sum = [sum(d2[t0:j+1]) for j in range(t0,len(d2))]
            ax.plot(t[i], d2)
            ax.plot(t[i], d1)
            ax.plot(t[i][t0:], d2sum)
            ax.fill_betweenx(acf[i], t[i][t0], t[i][t1], alpha=0.1)
            ax.axhline(2e-2)
        ax.set_yscale('log')
        ax.set_xlabel(r'$\delta t$')
        ax.set_ylabel(r'$ACF$')
        ax.set_ylim(1e-4, 1)
        ax.text(2*t[i][-1]/3, 0.5, r'$\tau=$'+str(int(tau)))
        ax.set_title(r'$\tilde{T}=$'+str(T))
    plt.show()
    return

def PlotFigTau():
    data = np.loadtxt('data/tauTemp1D.dat')
    plt.plot(data[:,0],data[:,1], '-o')
    plt.xlabel(r'$\tilde{T}$')
    plt.ylabel(r'$<\tau>$')
    plt.xlim(0,2000)
    plt.ylim(0, 2000)
    plt.show()

def main():
    # t, acf = LoadAcfData()
    # t, acf = ACFData()
    # PlotFig(t, acf)
    PlotFigTau()

if __name__ == "__main__":
    main()
        
