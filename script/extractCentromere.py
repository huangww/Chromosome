import numpy as np
import matplotlib.pyplot as plt
from theory import *

beadNumber = 18
ringSize = (beadNumber+2)/2
dimension = 3
cm = 5      # Centromere location
Teff = 1
dataDir = 'data/'

fileName = dataDir + 'MD_N' + str(beadNumber) + \
        '_T'+str(Teff) + '_sample.dat'
data = np.loadtxt(fileName,comments='#')
beadPos = data.reshape([-1,beadNumber,dimension])
distance = np.zeros((len(beadPos),ringSize,dimension))
for i in range(1,cm):
    distance[:,i,:]= beadPos[:,i+ringSize-1,:]-beadPos[:,i,:]
for i in range(cm,ringSize):
    distance[:,i,:]= beadPos[:,i+ringSize-2,:]-beadPos[:,i,:]
meanDis = distance.mean(axis=0)
varDis = distance.var(axis=0)

# save extracted data
fileName = fileName[:-4] + '_cm_xyz_varxyz.dat'
np.savetxt(fileName, np.vstack((meanDis.T,varDis.T)).T)

