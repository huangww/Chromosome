import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from theory import *
import glob
import os
import sys

figDir = 'fig/'
dataDir = sys.argv[1] + '/'
print dataDir
# dataDir = 'data/'

def GetParameter():
    fname = max(glob.iglob(dataDir + 'r*.dat'), key=os.path.getctime)
    parameters = fname.split('_')
    Nstr = [s[1:] for s in parameters if s[0] == 'N']
    Tstr = [s[1:] for s in parameters if s[0] == 'T']
    N = int(Nstr[0])
    T = float(Tstr[0])
    return N, T

def LoadData(N):
    dataFiles = glob.glob(dataDir + 'r*.dat')
    posList = []
    for dfile in dataFiles:
        data = np.loadtxt(dfile,comments='#')
        pos = data.reshape([-1,N,3])
        posList.append(pos)
    return np.array(posList)


def LoadMeanVar(N, T):
    fname = dataDir+'mean_N'+str(nBead)+'_T'+ '%g'%T +'.dat'
    data = np.loadtxt(fname)
    posMean = data.reshape([-1, N, 3])
    fname = dataDir+'var_N'+str(nBead)+'_T'+ '%g'%T +'.dat'
    data = np.loadtxt(fname)
    posVar = data.reshape([-1, N, 3])
    return posMean, posVar

def Save3DArray(fname, data):
    with file(fname, 'w') as f:
        f.write('# Data shape:{0}\n'.format(data.shape))
        t = 0
        for frame in data:
            f.write('# frame = '+ str(t) + '\n')
            np.savetxt(f, frame)
            t = t + 1

def UpdateFig(i, data1, l1, data2, l2):
    l1.set_ydata(data1[i,:,0])
    l2.set_ydata(data2[i,:,0])
    return l1, l2,

def MakeMovie(nBead, T, posMean, posVar):
    # mean and var from theory
    meanz = mean_z(nBead, T)
    varz = var_z(nBead, T)

    fig = plt.figure(0, figsize=(10,4))
    font = {'family' : 'scans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 12}
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    plt.rc('text', usetex = True)

    # generate a movie
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(meanz,'r-')
    ax2.plot(varz,'r-')
    l1, = ax1.plot(np.zeros(nBead), 'b-')
    l2, = ax2.plot(np.zeros(nBead), 'b-')

    ax1.set_xlabel('Bead index')
    ax1.set_ylabel(r"$<r_z>$")
    # ax1.set_ylim([-5,50])
    ax2.set_xlabel('Bead index')
    ax2.set_ylabel(r"$var[r_z]$")
    anima = animation.FuncAnimation(fig, UpdateFig, len(posMean), fargs=(posMean, l1, posVar, l2), interval=50)
    fname = figDir+'animaN'+str(nBead)+'T'+ '%g'%T +'.mp4'
    anima.save(fname)


def main(data = 'load'):
    nBead, T = GetParameter()
    # nBead, T = 100, 0.1
    if data == 'load':
        # load data and extract mean and var
        pos = LoadData(nBead)
        posMean = pos.mean(axis=0)
        posVar = pos.var(axis=0)

        # save mean and var
        fname = dataDir+'mean_N'+str(nBead)+'_T'+ '%g'%T +'.dat'
        Save3DArray(fname, posMean)
        fname = dataDir+'var_N'+str(nBead)+'_T'+ '%g'%T +'.dat'
        Save3DArray(fname, posVar)
    else:
        # Or, load mean and var directly
        posMean, posVar = LoadMeanVar(nBead, T)
    MakeMovie(nBead, T, posMean, posVar)

if __name__ == '__main__':
    main('load')
    # main('loaded')
    print 'done'
