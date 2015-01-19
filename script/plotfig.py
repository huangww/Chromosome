import matplotlib.pyplot as plt
import numpy as np
from theory import *

# initialize the figure
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 12 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text',usetex = True)
fig = plt.figure(0)
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

beadNumber = 300
# frameNumber = 5000
dimension = 3
index = np.linspace(0,beadNumber,beadNumber)
tempFrames = 0
step = 1

# for T in ['1', '5','10','5000']:
for T in [1,10,20,50,100]:
    # load data
    dataDir = 'data/'
    # fileName = dataDir + 'r_N' + str(beadNumber) + \
            # '_T'+str(T) + '_5489.dat'
    fileName = dataDir + 'MD_N' + str(beadNumber) + \
            '_T'+str(T) + '_sample.dat'
    data = np.loadtxt(fileName,comments='#')
    beadPos = data.reshape([-1,beadNumber,3])
    
    # calculate mean and varriance from data
    posMean = beadPos.mean(axis = 0)
    posVars = beadPos.var(axis = 0)

    # plot figure
    N = beadNumber
    line, = ax1.plot(index, z_mean(index,N,int(T)))
    ax1.plot(index, posMean[:,0],line.get_color()+'o')
    line, = ax2.plot(index, z_var(index,N,int(T)))
    ax2.plot(posVars[:,0],line.get_color()+'o')
    ax3.plot(index, posMean[:,1],line.get_color()+'o')
    ax3.plot(index, posMean[:,2],line.get_color()+'*')
    var_xy = np.zeros(len(index))
    for i in range(N):
        var_xy[i] = xy_var(i,N,int(T))
    line, = ax4.plot(index, var_xy/2.0)
    ax4.plot(index, posVars[:,1],line.get_color()+'o')
    ax4.plot(index, posVars[:,2],line.get_color()+'*')

ax3.set_ylim([-5,5])
ax1.set_xlabel('Bead index')
ax2.set_xlabel('Bead index')
ax3.set_xlabel('Bead index')
ax4.set_xlabel('Bead index')
ax1.set_ylabel(r"$<r_z>$")
ax2.set_ylabel(r"$var[r_z]$")
ax3.set_ylabel(r"$<r_{x,y}>$")
ax4.set_ylabel(r"$var[r_{x,y}]$")
plt.show()

