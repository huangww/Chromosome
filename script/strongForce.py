import numpy as np
import matplotlib.pyplot as plt
import theory as th
import acfFit as af

ax = plt.figure().add_subplot(111)

N = 100
Teff = np.linspace(0.2, 9.9, 98)

for T in Teff:
    fname = 'data/rg1D_N100_T%g_0.dat'%T
    data = np.loadtxt(fname)
    t = data[:,0]
    x = data[:,1]
    tau = af.GetTau(t, x)
    ax.scatter(T, tau)
    ax.scatter(T, T*(x.mean()))
# data = np.loadtxt('data/tau_N100_T10l.dat')
# ax.plot(Teff, data, 'o')

# zmeanT = [th.mean_z(N, T).mean() for T in Teff]
# tau = Teff * zmeanT
# zMid = [th.mean_z(N, T)[N/2] for T in Teff]
# tau = Teff * zMid
# ax.plot(Teff, tau)
ax.set_xlabel('1/F')
ax.set_ylabel('tau')
plt.show()
