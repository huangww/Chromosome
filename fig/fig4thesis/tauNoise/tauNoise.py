import numpy as np
import matplotlib.pyplot as plt

def tauTheory(F, T, L):
    """ calculate the theoretical relaxation time

    :F: strength of external force field
    :T: effective temperature
    :L: the total number of lattice sites
    :returns: the relaxation time tau

    """
    a, kb, xi = 1, 1, 1
    factor = np.exp((F*a)/(kb*T))
    tau = xi*a**2*L**2/(3*np.pi**2*kb*T) - xi*a**2*L**4*(factor-1)**2/(3*np.pi**2*kb*T*(L**2*(factor-1)**2+np.pi**2*(factor+1)**2))
    return tau

def plotTauNoise():
    """ plot the figure of relaxation time vs the noise level
    :returns: Null

    """
    fig = plt.figure(0,figsize=(5, 3.6))
    plt.rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 14 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.15, right =0.9,\
        bottom=0.15, top =0.9, wspace=0.25)

    # tau0 = (L/np.pi)**2/3
    T = np.linspace(1, 150, 100)
    ax.plot(T, tauTheory(1, T, 100), 'b-')
    ax.plot(T, tauTheory(0.1, T, 100), 'r-')
    # ax.plot(T, tauTheory(0.1, T, 100), '-')
    # ax.plot(T, tauTheory(0.01, T, 100), '-')

    # load the simulation data and plot as dots
    tauData = np.loadtxt('tauNoise_N100_F1.dat')
    ax.plot(tauData[:,0], tauData[:,1], 'bo')
    tauData = np.loadtxt('tauNoise_N100_F0.1.dat')
    ax.plot(tauData[:,0], tauData[:,1], 'ro')


    # set format of the figure
    ax.set_yscale('log')
    # ax.set_xlabel(r'$k_B T/\Delta E$')
    # ax.set_ylabel(r'$\tau_{\parallel}/\tau_{Rouse}$')
    ax.set_xlabel(r'$T$')
    ax.set_ylabel(r'$\tau_{\parallel}$')
    ax.legend([r'$F=1$', r'$F=0.1$'], loc='upper right', fontsize=15)
    ax.set_xlim([0,150])
    # ax.set_ylim([0,1.2])
    plt.show()
   


if __name__ == "__main__":
    # main()
    plotTauNoise()


    

