import glob as glob
import numpy as np
import matplotlib.pyplot as plt

def dataMerge(fname):
    fileList = glob.glob(fname)
    tauData = []
    for f in fileList:
        data = np.loadtxt(f)
        tauData.append(data)

    tau = np.nanmean(tauData, axis=0)
    print np.shape(tauData)
    print tau
    fname = '_'.join(fileList[0].split('_')[:-1])+'.dat'
    print fname
    # np.savetxt(fname, tau)
    
def PlotFig():
    fig = plt.figure(0,figsize=(5,4))
    plt.rc('text',usetex = True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 20 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    fig.subplots_adjust(left=0.2, right =0.95,\
        bottom=0.18, top =0.95, wspace=0.25)
    ax = fig.add_subplot(111)

    data1 = np.loadtxt('tau_N100_nPar10.dat')
    data2 = np.loadtxt('tau_N100_nPar20.dat')
    data3 = np.loadtxt('tau_N100_nPar50.dat')

    # T = np.linspace(1, 1000, 1000)
    # factor = np.exp(1/T)
    # p = np.sqrt(factor/(1+factor))
    # q = np.sqrt(1/(1+factor))
    # L = 100
    # tau = 1./((p-q)**2 + p*q*np.pi*np.pi/L**2)

    # ax.plot(T, tau)
    ax.plot(1./data1[:,0], data1[:,1]/2., 'o',label=r'$N=10$')
    ax.plot(1./data2[:,0], data2[:,1]/2., 's',label=r'$N=20$')
    ax.plot(1./data3[:,0], data3[:,1]/2., '*',label=r'$N=50$')
    # ax.legend(['theory', '# 10', '# 20', '# 50'], loc='lower right')
    ax.legend(loc='lower left',fontsize=20)
    ax.set_xscale('log')
    # ax.set_xlim([0.9,1.0])
    ax.set_xlabel(r'$q$',fontsize=20)
    # ax.set_xlabel(r'$\tilde{T}$')
    ax.set_ylabel(r'$\tau$',fontsize=20)
    # plt.savefig('diffnPar.pdf')

    plt.show()

def main():
    # fname = 'beadRodN100/*tau*T100_*.dat'
    # dataMerge(fname)
    PlotFig()

if __name__ == "__main__":
    main()
