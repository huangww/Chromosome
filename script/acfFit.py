import statsmodels.tsa.stattools as ss
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep, splev


dataDir = 'job.sh/'

def ACF(data):
    x = data[:, 1]
    t = data[:, 0]
    acf= ss.acf(x, nlags=len(x), fft=True)
    end = next(i for i,v in enumerate(acf) if v<0)
    return t[:end], acf[:end]

def GetACF(T):
    tlist = []
    acflist = []
    for i in range(1, 101):
        fname = dataDir + 'rg1D_N100_T'+str(T)+ \
                '_'+str(i)+'.dat'
        data = np.loadtxt(fname)
        t = data[:,0]
        x = data[:,1]
        acf = ss.acf(x, nlags=1e5, fft=True)
        tlist.append(t[:len(acf)])
        acflist.append(acf)
    acf = np.array(acflist).mean(axis=0)
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
    t13 = next((i for i,v in enumerate(d2sum) if v>2e-1), 
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
    ax.fill_betweenx(acf, t[t0], t[t1], alpha=0.1)
    tFit = t[t0:t1]
    acfFit = acf[t0:t1]
    coefs = np.polyfit(tFit, np.log(acfFit), 1)
    tau = -1./coefs[0]
    print 'tau =', tau
    ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]))
    ax.set_yscale('log')
    ax.set_xlabel(r'$\delta t$')
    ax.set_ylabel(r'$ACF$')
    ax.set_ylim(1e-3, 1)
    f = splrep(np.arange(len(acf)), np.log(acf))
    d1 = splev(np.arange(len(acf)), f, der=1)
    d2 = splev(np.arange(len(acf)), f, der=2)
    d2 = np.fabs(d2)
    d2sum = [sum(d2[t0:i+1]) for i in range(t0,len(d2))]
    ax.plot(t, d2)
    ax.plot(t, d1)
    ax.plot(t[t0:], d2sum)
    ax.axhline(3e-3)
    ax.axhline(2e-1)
    plt.show()
    fig.savefig('fig/acfFit.pdf')

def AddToFile(T, tau):
    f = open('data/tauTemp.dat', 'ab')
    np.savetxt(f, (T, tau), newline=' ')
    f.write('\n')
    f.close()

def LoadData(T):
    fname = dataDir + 'rg_N100_T'+str(T)+'_'+str(T)+'.dat'
    data = np.loadtxt(fname)
    return data

def main():
    T = 20
    t, acf = GetACF(T)
    fname = 'data/acf_N100_T'+str(T)+'.dat'
    np.savetxt(fname, (t, acf))
    # t, acf = np.loadtxt(fname)
    tau = FitTau(t, acf)
    PlotFig(t, acf)
    AddToFile(T, tau)
     
if __name__ == "__main__":
    main()
