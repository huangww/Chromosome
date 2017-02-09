import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

def qNumber(x, q):
    if q == 1:
        return x
    else:
        return (1-q**x)/(1-q)

def qBracket(x, q):
    return (q**x - q**(-x) ) / ( q - 1./q)

def qBinorm(N, x, q):
    numerator = 1
    denominator = 1
    for i in range(1, x+1):
        # denominator = denominator * qBracket(i, q)
        # numerator = numerator * qBracket(N-i+1, q)
        denominator = denominator * (1-q**(i))
        numerator = numerator * (1-q**(N-i+1))
    return numerator / denominator

def tagPdf(T, k, N, L):
    factor = np.exp(-1./T)
    p = factor/(1+factor)
    q = 1./(1+factor)
    x = np.arange(1, L+1)
    q = p / q
    # q = 1
    qx = np.zeros(L)
    for i in range(L):
        qx[i] = qBinorm(x[i]-1, k-1, q) * qBinorm(L-x[i], N-k, q) 
    px = q**((x-k)*(N+1-k)) * qx / qBinorm(L, N, q)
    return px

def density(T, N, L):
    rho = np.zeros(L)
    for k in range(1, N+1):
        rho = rho + tagPdf(T, k, N, L)
    return rho

def densitySchutz(T, N, L):
    factor = np.exp(-1./T)
    p = factor/(1+factor)
    q = 1./(1+factor)
    q = np.sqrt(p/q)
    x = np.arange(1, L+1)
    rho = []
    rho = q**(-L-1+2*x) / qBracket(L, q)
    for n in range(2, N+1):
        rho = qBracket(n, q) * q**(-L-1+2*x) * abs(1 - rho) / qBracket(L-n+1, q)
    return np.array(rho)


def plotCompareFig():
    plt.close('all')
    plt.rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 15 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    fig = plt.figure(0, figsize=(6, 4.2))
    fig.subplots_adjust(left=0.13, right =0.95,\
            bottom=0.13, top =0.95, wspace=0.25)
    ax = fig.add_subplot(111)
    # px = densitySchutz(T, N, L)
    # px = tagPdf(T, 30, N, L)
    Tarray = [1, 10, 100, 1000]
    N, L = (50, 100)
    x = np.arange(1, L+1)
    for T in Tarray:
        px = density(T, N, L)
        fx = 1./(1.+np.exp((x-(L+1)/2.)/float(T)))
        labelString = r"$\frac{k_B T}{\Delta E}=$" + \
                r"${}$".format(T)
        line, = ax.plot(x, px, '-', label=labelString)
        ax.plot(x, fx, 'o', c = line.get_color())
    ax.legend(fontsize=15)
    ax.set_xlabel(r'$\rm{Index}$')
    ax.set_ylabel(r'$\mathbb{P}(n_j)$')
    axinset = inset_axes(ax, width="40%", height="40%", loc=3)
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    Tarray = [1, 10, 100]
    N, L = (5, 10)
    x = np.arange(1, L+1)
    for T in Tarray:
        px = density(T, N, L)
        fx = 1./(1.+np.exp((x-(L+1)/2.)/float(T)))
        line, = axinset.plot(x, px, '-', label='1/F='+str(T))
        axinset.plot(x, fx, 'o', c = line.get_color())
    plt.show()

def plotCompareFig2():
    plt.close('all')
    plt.rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 15 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    fig = plt.figure(0, figsize=(10, 5))
    Tarray = [1, 2, 5, 10, 1000]
    labelStr0 = r"$\tilde{T}=$"
    labelStr = [r"$1$", r"$2$", r"$5$", r"$10$", r"$10^3$"]
    labelStr = [labelStr0 + s for s in labelStr]
    N, L = (5, 10)
    x = np.arange(1, L+1)
    ax = fig.add_subplot(121)
    for i,T in enumerate(Tarray):
        px = density(T, N, L)
        fx = 1./(1.+np.exp((x-(L+1)/2.)/float(T)))
        line, = ax.plot(x, px, '-', label=labelStr[i])
        ax.plot(x, fx, 'o', c = line.get_color())
    ax.set_xlabel(r'$\rm{Index}$')
    ax.set_ylabel(r'$\mathbb{P}(n_j)$')
    ax.legend(loc='lower left',fontsize=15)

    Tarray = [1, 10, 50, 100, 1000]
    labelStr0 = r"$\tilde{T}=$"
    labelStr = [r"$1$", r"$10$", r"$50$", r"$10^2$", r"$10^3$"]
    labelStr = [labelStr0 + s for s in labelStr]
    N, L = (50, 100)
    x = np.arange(1, L+1)
    ax = fig.add_subplot(122)
    for i,T in enumerate(Tarray):
        px = density(T, N, L)
        fx = 1./(1.+np.exp((x-(L+1)/2.)/float(T)))
        labelString = r"$\tilde{T}$=" + \
                r"${}$".format(str())
        line, = ax.plot(x, px, '-', label=labelStr[i])
        ax.plot(x, fx, 'o', c = line.get_color())
    ax.set_xlabel(r'$\rm{Index}$')
    ax.set_ylabel(r'$\mathbb{P}(n_j)$')
    ax.legend(loc='lower left', fontsize=15)
    plt.show()

def main():
    plotCompareFig2()
    return 

if __name__ == "__main__":
    main()
