import numpy as np
import matplotlib.pyplot as plt
from theory import *

beadNumber = 598
ringSize = (beadNumber+2)/2
# ringSize = (beadNumber+1)/2
dimension = 3
cm =  150     # Centromere location
Teff = 100
dataDir = 'data/'

fileName = dataDir + 'MD_N' + str(beadNumber) + \
        '_T'+str(Teff) + '_sample.dat'
        # '_T'+str(Teff) + '.dat'
data = np.loadtxt(fileName,comments='#')
beadPos = data.reshape([-1,beadNumber,dimension])
distance = np.zeros((len(beadPos),ringSize+1,dimension))

# centromere rigid bond
for i in range(1,cm):
    distance[:,i,:]= beadPos[:,i+ringSize-1,:]-beadPos[:,i,:]
for i in range(cm,ringSize):
    distance[:,i,:]= beadPos[:,i+ringSize-2,:]-beadPos[:,i,:]

# centromere bond with spring
# for i in range(1,ringSize):
#     distance[:,i,:]= beadPos[:,i+ringSize-1,:]-beadPos[:,i,:]

meanDis = distance.mean(axis=0)
varDis = distance.var(axis=0)

# save extracted data
fileName = fileName[:-4] + '_cm_xyz_varxyz.dat'
np.savetxt(fileName, np.vstack((meanDis.T,varDis.T)).T)

