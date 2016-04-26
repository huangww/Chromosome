import numpy as np
import matplotlib.pyplot as plt


def FermiDrac(T, N):
    mu = (N-1)/2.
    x = np.arange(N)
    p = 1. / (1. + np.exp((x-mu)/T))
    x = 2.*x/N
    return x, p
    
    
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
        
def main():
    Demo()

if __name__ == "__main__":
    main()
