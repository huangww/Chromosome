import numpy as np
import os


beadNumber = 500
dimension = 3
Teff = 1

sampleStep = 10

dataDir = 'data/'
fileName = dataDir + 'MD_N' + str(beadNumber) + '_T' + str(Teff) + '_sample.dat'
sampleFile = dataDir + 'MD_N' + str(beadNumber) + '_T' + str(Teff) + '_sample.dat'


# solution 1:

beadNumber = beadNumber + 1
sample = open(sampleFile, 'w')
count = 0
with open(fileName) as f:
    for line in f:
        if count % (beadNumber * sampleStep) < beadNumber:
            sample.write(line)
        count = count + 1

sample.close()
# os.rename(sampleFile, dataFile)



# solution 2:

# data = np.loadtxt(fileName, comments = '#')
# lineNumber = len(data[:,0])
# sampleData = np.zeros([lineNumber/sampleStep, dimension])
# count = 0
# for line in range(lineNumber):
    # if line % (beadNumber * sampleStep) < beadNumber:
        # sampleData[count] = data[line]
        # count = count + 1

# np.savetxt(fileName, sampleData, fmt = '%f', delimiter = '\t')