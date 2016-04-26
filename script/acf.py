import numpy as np
import matplotlib.pyplot as plt
import statsmodels.tsa.stattools as ss

dataDir = 'data/'

def LoadData(N, T):
    fname = dataDir + 'r_N'+str(N)+'_T'+str(T)+'_0.dat'
    data = np.loadtxt(fname)
    pos = data.reshape([-1,N,3])
    return pos

def ACF(x):
    xvar = x.var()
    x = x - x.mean()
    n = len(x)
    lag = np.arange(100)
    acf = np.zeros(len(lag))
    acf = np.array([(x[:n-i] * x[i:]).sum() / (xvar*(n-i))
        for i in lag])
    return acf

def main():
    N = 100
    T = 1
    pos = LoadData(N, T)
    x = pos[:, 50, 0]
    acf = ACF(x)
    t = np.arange(len(acf))*1
    plt.plot(t, acf)
    plt.plot(t, ss.acf(x, nlags=99))
    plt.show()

if __name__ == "__main__":
    main()

