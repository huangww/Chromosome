import os as os


dataDir = './data/lammps/N50T1LJ/'

for files in os.listdir(dataDir):
    if files[:4] == 'dump':
        fname = 'dump' + str(int(files[4:-4])/100000) + '.dat'
        os.rename(dataDir+files, dataDir + fname)

