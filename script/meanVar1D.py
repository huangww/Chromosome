import numpy as np
import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from theory import *

N = 100
T = 1

figDir = 'fig/'
dataDir = 'data/'

def LoadData(T):
    dataDir = 'data/'
    dfile = 'par_N'+str(N)+'_T'+str(T)+'.dat'
    print dataDir + dfile
    data = np.loadtxt(dataDir + dfile)
    return np.array(data[:,1:]), np.array(data[:,0])
    
def PosMeanVar(pos0, t0):
    nbins = 1000
    binrange = (0,1000)
    posMean = np.zeros([nbins, N])
    posVar = np.zeros([nbins, N])
    binCount, t = np.histogram(t0, bins=nbins, 
            range=binrange)
    for i in range(len(pos0[0])):
        hist1, t = np.histogram(t0, bins=nbins, 
                range=binrange, weights=pos0[:,i])
        posMean[:,i]  = hist1 / binCount
        hist2, t = np.histogram(t0, bins=nbins, 
                range=binrange, weights=pos0[:,i]*pos0[:,i])
        posVar[:,i] = hist2 / binCount - posMean[:,i]*posMean[:,i]
    return posMean, posVar, t[:-1]

def PosMeanVarMD(pos0, t0):
    nbins = 1000
    pos = pos0.reshape([-1,nbins,len(pos0[0])])
    posMean = pos.mean(axis=0)
    posVar = pos.var(axis=0)
    t = t0[:nbins]
    return posMean, posVar, t

def UpdateFig(i, data1, l1, data2, l2):
    l1.set_ydata(data1[i])
    l2.set_ydata(data2[i])
    return l1, l2,

def MakeMovie(N, T, posMean, posVar):
    meanz = mean_z_1D_sum(N,2*T)
    varz = var_z_1D(N,2*T)

    fig = plt.figure(0, figsize=(10,4))
    font = {'family' : 'scans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 12}
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    plt.rc('text', usetex = True)

    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(meanz,'r-')
    ax2.plot(varz,'r-')
    l1, = ax1.plot(posMean[0], 'b-')
    l2, = ax2.plot(posVar[0], 'b-')

    ax1.set_xlabel('Bead index')
    ax1.set_ylabel(r"$<r_z>$")
    ax2.set_xlabel('Bead index')
    ax2.set_ylabel(r"$var[r_z]$")
    anima = animation.FuncAnimation(fig, UpdateFig, len(posMean), fargs=(posMean, l1, posVar, l2), interval=50)
    anima.save(figDir + 'anima1DN'+str(N) +'T'+str(T)+'.mp4')

def main():
    pos0, t0 = LoadData(T)
    posMean, posVar, t = PosMeanVar(pos0, t0)
    # posMean, posVar, t = PosMeanVarMD(pos0, t0)

    # save mean and var data
    fname = dataDir + 'asep_mean_N' + str(N) \
           + '_T' + str(T) + '.dat'
    np.savetxt(fname, np.vstack((t.T, posMean.T)).T)
    fname = dataDir + 'asep_var_N' + str(N) \
             + '_T' + str(T) + '.dat'
    np.savetxt(fname, np.vstack((t.T, posVar.T)).T)
    MakeMovie(N, T, posMean, posVar)


if __name__ == '__main__':
    main()
    print 'done'
