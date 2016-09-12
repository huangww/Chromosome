import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp


def FermiDrac(T, N):
    mu = (N+1)/2.
    x = np.arange(N)
    p = 1. / (1. + np.exp((x-mu)/T))
    x = 2.*x/N
    return x, p
    
def DensityLattice(T, N):
    a, F, t0, k = (1, 1, 1, 1)
    L = 2*N
    alpha = 1/(t0*(1+np.exp(-2*F*a/(k*T))))
    beta = np.exp(-2*F*a/(k*T))/(t0*(1+np.exp(-2*F*a/(k*T))))
    q = np.sqrt(alpha/beta)
    qfacNorm = mp.qfac(L, q)/(mp.qfac(N, q)*mp.qfac(N, q))
    rho = np.zeros(N)
    for x in range(1, N+1):
        sumTerm = 0
        for k in range(N):
            qfacK = mp.qfac(L, q)/(mp.qfac(L-x, q)*mp.qfac(x, q))
            sumTerm = sumTerm + (-1)**(N-k+1)*q**(-(N-k)*(L+1-2*x))*qfacK
        rho[x-1] = sumTerm
    rho = rho / qfacNorm
    return rho


def Demo():
    plt.close('all')
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    Tarray = [500, 250, 50, 5]
    N = 1000
    for T in Tarray:
        x, p = FermiDrac(T, N)
        ax.plot(x, p)
    ax.set_xlabel(r'$\epsilon_i/\mu$', fontsize=20, labelpad=1)
    ax.set_ylabel(r'$<n_i>$', fontsize=20, labelpad=1)
    ax.ticklabel_format(fontsize=20)
    ax.legend(['$k_B T = \mu$', '$k_B T = \mu/2$',\
            '$k_B T = \mu/10$', '$k_B T = \mu/100$'])
    plt.show()
        
def DemoLattice():
    plt.close('all')
    fig = plt.figure(0)
    ax = fig.add_subplot(111)
    ax.plot(DensityLattice(1000., 100))
    ax.set_xlabel('index')
    ax.set_ylabel(r'$<n_i>$')
    plt.show()

def main():
    # Demo()
    DemoLattice()

if __name__ == "__main__":
    main()
