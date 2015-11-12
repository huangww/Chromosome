import numpy as np
import matplotlib.pyplot as plt
from theoryGyration import *

font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 16 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text',usetex = True)

fig = plt.figure(0,figsize=(5,4))
fig.subplots_adjust(left=0.15, right =0.95,\
        bottom=0.15, top =0.95, wspace=0.25)
ax = fig.add_subplot(111)

dataDir = 'data/rg/'

beadNumber = 100
for T in [1, 2, 5, 10, 20, 50, \
        100, 200, 500, 1000, 2000, 5000]: 
    fileName = dataDir + 'rg_N' + str(beadNumber) + \
            '_T'+str(T) + '_5489.dat'
    data = np.loadtxt(fileName,comments='#')
    # ax.plot(data)

    gyration = data.mean()
    ax.scatter(T, gyration,\
            marker = 'o', color = 'g')

beadNumber = 100
Tlog = np.logspace(0, 4, 100)
Rgs = np.zeros(len(Tlog))
count = 0
for T in Tlog:
    Rgs[count] = Gyration(beadNumber, T)
    count = count + 1

line0, = ax.plot(Tlog, Rgs, 'b-')
line1, = ax.plot([1,10000], [beadNumber/12.0, beadNumber/12.0], 'r-')
# ax.legend((line0,line1),('Theory','Random Ring'))
ax.legend((line0, line1),('Theory', 'High $T$ asymptotic'))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([1,1e4])
ax.set_ylim([5,2e2])
ax.set_xlabel(r'T')
ax.set_ylabel(r"$\left<R_g^2\right>$")
plt.show()
