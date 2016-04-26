import numpy as np
import matplotlib.pyplot as plt
import acfFit as af



# Nbead = [10, 20, 30, 50, 70, 100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000, 10000]
Nbead = [10, 20, 30, 50, 70, 100, 200, 300, 500]

ax = plt.figure(0).add_subplot(111)

tau = []
for N in Nbead:
    fname = "data/rg1D_N" + str(N) + "_T2000_0.dat"
    data = np.loadtxt(fname)
    t = data[:,0]
    x = data[:,1]
    tauVal = af.GetTau(t, x)
    tau.append(tauVal)
ax.plot(Nbead, tau, 'o')

NbeadFit = np.log(Nbead)
tauFit = np.log(tau)
coefs = np.polyfit(NbeadFit, tauFit, 1)
exponent = coefs[0]

N = np.logspace(1, 3, 100)
ax.plot(N, np.exp(coefs[0]*np.log(N)+coefs[1]))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('N')
ax.set_ylabel('tau')
ax.set_title(str(exponent))

plt.show()
