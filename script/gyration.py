import numpy as np
import matplotlib.pyplot as plt
from theory import *

font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 12 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text',usetex = True)

fig = plt.figure(0)
ax = fig.add_subplot(111)

beadNumber = 500

dataDir = 'data/rg/'
for T in [1, 2, 5, 10, 20, 50, \
        100, 200, 500, 1000, 2000, 5000]: 
    fileName = dataDir + 'rg_N' + str(beadNumber) + \
            '_T'+str(T) + '_5489.dat'
    data = np.loadtxt(fileName,comments='#')
    # ax.plot(data)

    gyration = data.mean()
    ax.scatter(T, gyration,\
            marker = 'o', color = 'g')

Tlog = np.logspace(0, 4, 100)
Rgs = np.zeros(len(Tlog))
count = 0
for T in Tlog:
    Rgs[count] = Gyration(beadNumber, T)
    count = count + 1

ax.plot(Tlog, Rgs, 'b-')
ax.plot([1,10000], [beadNumber/12.0, beadNumber/12.0], 'r-')

ax.set_xscale('log')
ax.set_yscale('log')
# ax.set_xlim([1,10000])
# ax.set_ylim([10,10000])
# ax.set_xlabel('Time')
ax.set_xlabel('Dimensionless temperature')
ax.set_ylabel(r"$R_g^2$")
plt.show()
