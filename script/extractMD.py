import numpy as np
import matplotlib.pyplot as plt


dataDir = 'data/'
beadNumber = 300
dimension = 3
Teff = [1,10,20,50,100]

for T in Teff:
    fileName = dataDir + 'MD_N' + str(beadNumber) + '_T' + str(T) + '_sample.dat'
    data = np.loadtxt(fileName, comments= '#')

    beadPos = data.reshape([-1, beadNumber, dimension])


    meanPos = beadPos.mean(axis=0)
    varPos = beadPos.var(axis=0)


    plt.plot(meanPos[:,0])
    fileName = fileName[:-4] + '_xyz_varxyz.dat'

    np.savetxt(fileName, np.vstack((meanPos.T,varPos.T)).T)
