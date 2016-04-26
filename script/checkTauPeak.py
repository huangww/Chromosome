import numpy as np
import matplotlib.pyplot as plt
import acfFit as af



N = 100
Teff = np.concatenate([np.linspace(1, 9.5, 18), np.linspace(10, 95, 18), np.linspace(100, 1000, 19)])
# Teff = [10, 50, 100, 500, 1000]
# Teff = np.linspace(1, 10, 19)
ax = plt.figure().add_subplot(111)

for T in Teff:
    fname = "data/rg1D_N100_T%g_0.dat"%T
    data = np.loadtxt(fname)
    t = data[:,0]
    x = data[:,3]
    print T
    tau = af.GetTau(t, x)
    ax.scatter(T, tau, color='r')
    fname = "data/par_N100_T%g_0.dat"%T
    data = np.loadtxt(fname)
    pos = data[:,0]
    tau = af.GetTau(t, pos)
    ax.scatter(T, tau, color='b')
    pos = data[:,1]
    tau = af.GetTau(t, pos)
    ax.scatter(T, tau, color='g')
    # acf = af.ACF(x)
    # t = t[:len(acf)]
    # tau = af.FitTau(t, acf)
    # fname = "fig/acf_N100_T%g.pdf"%T
    # af.PlotFig(t, acf, fname)

ax.set_xlabel('1/F')
ax.set_ylabel('tau')

plt.show()
