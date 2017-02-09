import numpy as np
import matplotlib.pyplot as plt

def tauTheory(L):
    tau0 = (L/np.pi)**2
    Teff = np.linspace(1, 1000, 1000)
    factor = np.exp(-1./Teff)
    a = 2.*factor/(1.+factor)
    b = 2./(1.+factor)
    tau = 1. / (a + b - 2.*np.sqrt(a*b)*np.cos(np.pi/L))
    return tau/tau0

def tauAtT(T, L):
    # tau0 = (L/np.pi)**2
    T = float(T)
    factor = np.exp(-1./T)
    a = 2.*factor/(1.+factor)
    b = 2./(1.+factor)
    tau = 1. / (a + b - 2.*np.sqrt(a*b)*np.cos(np.pi/L))
    return tau

def plotTauLength():
    fig = plt.figure(0)
    plt.rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 14 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    ax = fig.add_subplot(111)

    Teff = [10, 20, 100, 10000]
    length = np.linspace(10, 10000, 1000)
    for T in Teff:
        tau = [tauAtT(T, l) for l in length]
        ax.plot(length, tau, label='F='+str(1./T))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Size')
    ax.set_ylabel(r'$\tau$')
    # ax.set_xlim([0,1000])
    # ax.set_ylim([0,1.2])
    ax.legend(loc='upper left', fontsize=14)

    plt.show()


def plotTauForce():
    fig = plt.figure(0)
    plt.rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 14 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    ax = fig.add_subplot(111)

    length = [100, 200, 500, 1000]
    for l in length:
        tau = tauTheory(l)
        ax.plot(tau, label='L='+str(l))

    ax.set_xscale('log')
    ax.set_xlabel(r'$k_B T / \Delta E$')
    ax.set_ylabel(r'$\tau / \tau_{Rouse}$')
    ax.set_xlim([0,1000])
    ax.set_ylim([0,1.2])
    ax.legend(loc='upper left', fontsize=14)

    plt.show()

def main():
    # plotTauForce()
    plotTauLength()

if __name__ == "__main__":
    main()

