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

# simulation data
dataDir = 'data/rg/'
size = np.array([10, 20, 50, 100, 200, 500, 1000, 2000])
for T,m,c in zip([1, 50, 1000],['o','*','s'],['r','b','g']):
    for N in size: 
        fileName = dataDir + 'rg_N' + str(N) + \
                '_T'+str(T) + '_5489.dat'
        data = np.loadtxt(fileName,comments='#')

        gyration = data.mean()
        ax.scatter(int(N), gyration, marker = m, color = c)

# theory
    Nlog = np.logspace(1, 4, 20)
    Rgs = np.zeros(len(Nlog))
    count = 0
    for N in Nlog:
        Rgs[count] = Gyration(int(N/2)*2, T)
        print int(N/2)*2, T, Rgs[count]
        count = count + 1
        
    ax.plot(Nlog, Rgs, c+'-')
    np.savetxt('rgN'+str(T)+'.dat',np.array([Nlog, Rgs]).T)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([1,1e4])
ax.set_ylim([10,1e4])
ax.set_xlabel(r'N')
ax.set_ylabel(r"$R_g^2$")
plt.show()
