import numpy as np
import os


dataDir = 'data/'
fileName = dataDir + 'timming.dat'
outputFile = fileName[:-4] + '_ed.dat'


timeFile = open(outputFile, 'w')

with open(fileName) as f:
    for line in f:
        i0 = line.find('.')
        i0 = line.find('.',i0+1)
        i0 = line.find('.',i0+1)
        j0 = line.find(':')
        i = line.find(' ')
        j = line.find(' ', i+1)
        line = line[i0+1:j0]+'\t'+line[i+1:j]+'\n'
        timeFile.write(line)

timeFile.close()
