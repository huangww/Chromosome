import matplotlib.pyplot as plt
import numpy as np

beadNumber = 50
frameNumber = 10000
dimension = 3
index = np.linspace(0,beadNumber,beadNumber)
tempFrames = 10000
step = 1

# initialize the figure
fig = plt.figure(0)
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)

for Temp in ['1', '10']:
    # load data
    dataDir = 'data/lammps/N50T' + Temp + 'LJ/'
    beadPos = np.zeros((frameNumber, beadNumber, dimension))
    t = 0
    while t<frameNumber:
        datafile = dataDir + 'dump' + str((t+tempFrames)*step) + '.dat'
        posFrame = np.loadtxt(datafile, skiprows = 2, usecols = (1, 2, 3))
        beadPos[t] = posFrame
        t = t + 1

    # calculate mean and varriance from data
    posMean = beadPos.mean(axis = 0)
    posVars = beadPos.var(axis = 0)
 
    # load MD data
    dataDir = 'data/MD_N50/'
    # determine temperature from file names
    if Temp == '0':
        T = '50'
    else:
        T = Temp
    datafile = dataDir + 'MD_N50_T'+T+'_xyz_varxyz_pseudo_LJ.dat'
    dataMD = np.loadtxt(datafile)

    # plot figure
    line, = ax1.plot(index, dataMD[0,:])
    ax1.plot(index, posMean[:,0],line.get_color()+'o')
    line, = ax2.plot(index, dataMD[3,:])
    ax2.plot(posVars[:,0],line.get_color()+'o')
    line, = ax3.plot(index, (dataMD[1,:]+dataMD[2,:])/2)
    ax3.plot(index, posMean[:,1],line.get_color()+'o')
    ax3.plot(index, posMean[:,2],line.get_color()+'*')
    line, = ax4.plot(index, (dataMD[4,:]+dataMD[5,:])/2)
    ax4.plot(index, posVars[:,1],line.get_color()+'o')
    ax4.plot(index, posVars[:,1],line.get_color()+'*')

ax3.set_ylim([-5,5])

plt.show()

