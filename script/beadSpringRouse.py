import numpy as np
import matplotlib.pyplot as plt
import itertools as it
import acfFit as af

def RouseMean(N, F):
    kh = 3.
    prefactor = F/(2.*N*kh)
    z = np.zeros(N)
    for k in range(N):
        z[k] = 0
        sumterm = 0
        for (m,j) in it.product(range(1,N), range(1,N)):
            sumterm = sumterm + np.sin(k*m*np.pi/N)* \
            np.sin(j*m*np.pi/N)/(np.sin(m*np.pi/(2.*N))**2)
        z[k] = prefactor * sumterm
    return z


def RouseVar(N, F):
    kh = 3.
    prefactor = 3./(2*N*kh)
    var = np.zeros(N)
    for k in range(N):
        var[k] = 0
        sumterm = 0
        for m in range(1, N):
            sumterm = sumterm + (np.sin(k*m*np.pi/N)/ \
                    np.sin(m*np.pi/(2.*N)))**2
        var[k] = prefactor * sumterm
    return var
    
def RouseAcf(N, k, t):
    kh = 3
    prefactor = 3./(2.*N*kh)
    sumterm = 0
    for m in range(1, N):
        sumterm = sumterm + (np.sin(k*m*np.pi/N)/ \
                np.sin(m*np.pi/(2.*N)))**2 * \
                np.exp(-kh*4*np.sin(m*np.pi/(2.*N))**2*t) 
    acf = prefactor * sumterm
    acf = acf / acf.max()
    return acf



def PlotFigMean():
    N = 100
    ax = plt.figure().add_subplot(111)
    Fext = [1., 0.5, 0.2, 0.1, 0.05]
    # Fext = [1.]
    index = np.arange(N)
    for F in Fext:
        meanz = RouseMean(N, F)
        line, = ax.plot(index, meanz, lw=2, label='F='+str(F))
        fname = 'data/rMeanVar_N%g_T%g.dat'%(N, 1/F)
        data = np.loadtxt(fname)
        ax.plot(index[::5], data[0,::5], 'o', c=line.get_color())
    ax.set_xlabel('Index', fontsize=20)
    ax.set_ylabel(r'$z/a$', fontsize=20)
    ax.legend()
    plt.show()

def PlotFigVar():
    N = 100
    ax = plt.figure().add_subplot(111)
    Fext = [1., 0.5, 0.2, 0.1, 0.05]
    # Fext = [1.]
    index = np.arange(N)
    var = RouseVar(N, 0)
    for F in Fext:
        line, = ax.plot(index, var, lw=2, label='F='+str(F))
        fname = 'data/rMeanVar_N%g_T%g.dat'%(N, 1/F)
        data = np.loadtxt(fname)
        ax.plot(index[::5], 3*data[-1,::5], 'o', c=line.get_color())
    ax.set_xlabel(r'$\rm{Index}$', fontsize=20)
    ax.set_ylabel(r'$\rm{var}[\mathbf{r}]$',fontsize=20)
    ax.legend()
    plt.show()

def PlotAcf(tid):
    N = 100
    ax = plt.figure().add_subplot(111)
    index = [10, 20, 50]
    fname = 'rd_N100_T10_%g.dat'%tid
    pos = np.loadtxt(fname).reshape([-1,3,3])
    fname = fname.replace('rd_', 'rg_')
    t = np.loadtxt(fname)[:,0]
    for i in index:
        acf = af.ACF(pos[:,index.index(i),2])
        line, = ax.plot(t[:len(acf)], acf)
        t0 = np.linspace(0, t[len(acf)], 1000)
        acf = RouseAcf(N, i, t0)
        ax.plot(t0, acf, '--', lw=2, c=line.get_color(), \
                label=r'$i=$'+str(i))
    ax.set_xlabel(r'$t$', fontsize=20)
    ax.set_ylabel(r'$ACF$',fontsize=20)
    ax.set_yscale('log')
    ax.legend()
    fname = fname.replace('.dat', '.pdf')
    fname = fname.replace('rg_', 'acf_')
    print fname
    plt.savefig(fname)
    plt.show()

def PlotAcf1D():
    N = 100
    ax = plt.figure().add_subplot(111)
    fname = 'data/rg1D_N100_T100_1.dat'
    data = np.loadtxt(fname)
    t = data[:,0]
    for i in range(1, len(data[0])):
        acf = af.ACF(data[:,i])
        line, = ax.plot(t[:len(acf)], acf)
        t0, t1 = af.FindFitRange(acf)
        tFit = t[t0:t1]
        acfFit = acf[t0:t1]
        coefs = np.polyfit(tFit, np.log(acfFit), 1)
        ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]), \
                c=line.get_color(), lw=2)
    ax.set_xlabel(r'$t$', fontsize=20)
    ax.set_ylabel(r'$ACF$',fontsize=20)
    ax.set_yscale('log')
    ax.legend(['Rg', '10', '20', '30', '40', '50'])
    fname = fname.replace('.dat', '.pdf')
    fname = fname.replace('rg1D_', 'acf1D_')
    print fname
    plt.savefig(fname)
    plt.show()

def PlotAcfRg():
    N = 100
    ax = plt.figure().add_subplot(111)
    Teff = [50, 100,150, 200]
    for T in Teff:
        fname = 'data/rg1D_N100_T%g_1.dat'%T
        data = np.loadtxt(fname)
        t = data[:,0]
        rg = data[:,1]
        acf = af.ACF(rg)
        line, = ax.plot(t[:len(acf)], acf, label=str(T))
        t0, t1 = af.FindFitRange(acf)
        tFit = t[t0:t1]
        acfFit = acf[t0:t1]
        coefs = np.polyfit(tFit, np.log(acfFit), 1)
        ax.plot(tFit, np.exp(coefs[0]*tFit+coefs[1]), \
                c=line.get_color(), lw=2)
    ax.set_xlabel(r'$t$', fontsize=20)
    ax.set_ylabel(r'$ACF$',fontsize=20)
    ax.set_yscale('log')
    ax.legend()
    fname = 'acf_rg_N100.pdf'
    plt.savefig(fname)
    plt.show()


def main():
    # PlotFigMean()
    # PlotFigVar()
    # for i in range(1,101):
    #     PlotAcf(i)
    PlotAcf1D()
    # PlotAcfRg()
    

if __name__ == "__main__":
    main()
