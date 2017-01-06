import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import control
import fileinput
import time

def msd_straight_forward(r):
    shifts = np.arange(len(r))
    msds = np.zeros(shifts.size)    
    for i, shift in enumerate(shifts):
        diffs = r[:-shift if shift else None] - r[shift:]
        sqdist = np.square(diffs)
        msds[i] = sqdist.mean()
    return msds

def autocorrFFT(x):
    N=len(x)
    F = np.fft.fft(x, n=2*N)  #2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res= (res[:N]).real   #now we have the autocorrelation in convention B
    n=N*np.ones(N)-np.arange(0,N) #divide res(m) by (N-m)
    return res/n #this is the autocorrelation in convention A

def msd_fft(r):
    N=len(r)
    # D=np.square(r).sum(axis=1) 
    D=np.square(r)
    D=np.append(D,0) 
    # S2=sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    S2 = autocorrFFT(r)
    Q=2*D.sum()
    S1=np.zeros(N)
    for m in range(N):
      Q=Q-D[m-1]-D[N-m]
      S1[m]=Q/(N-m)
    return S1-2*S2

def getMSD(fname):
    # fname = max(glob.iglob('data/BeadRod_r_N2_*.dat'), \
    #         key=os.path.getctime)
    pos = np.loadtxt(fname).reshape([-1,2,2])
    vector = pos[:,1,:2] - pos[:,0,:2]
    vectorNorm = np.sqrt(np.sum(vector**2, axis=1))
    cosang = vector[:,0] / vectorNorm
    # sinang = np.sqrt(np.sum(vector[:,1:]**2, axis=1))
    sinang = vector[:,1] / vectorNorm
    theta = np.arctan2(sinang, cosang)
    fname = fname.replace('_r_', '_rd_')
    t = np.loadtxt(fname)[:,0]
    theta = control.unwrap(theta+np.pi)-np.pi
    msd = msd_fft(theta)
    return t, msd

def averagedMSD():
    fname = 'data/BeadRod_r_N2_*.dat'
    fileList = glob.glob(fname)
    t, msd = getMSD(fileList[0])
    for f in fileList[1:]:
        t, msd0 = getMSD(f)
        msd = msd + msd0
    msd = msd / len(fileList)    
    return t, msd

def plotFig():
    # t, msd = averagedMSD()
    fname = 'data/msd_2D_dumbbell.dat'
    # np.savetxt(fname, [t, msd])
    t, msd = np.loadtxt(fname)
    coefs = np.polyfit(t, msd, 1)
    plt.plot(t, coefs[0]*t+coefs[1],'r')
    plt.plot(t, msd, label='dumbbell 2D')
    print coefs
    # fname = 'data/msd_dumbbell.dat'
    fname = 'data/msd_2D_pinnedBead.dat'
    # np.savetxt(fname, [t, msd])
    t, msd = np.loadtxt(fname)
    coefs = np.polyfit(t, msd, 1)
    plt.plot(t, coefs[0]*t+coefs[1],'r')
    plt.plot(t, msd, label='pinned bead 2D')
    print coefs
    plt.xlabel(r'$t$')
    plt.ylabel(r'MSD$(\theta)$')
    plt.legend(loc='upper left')
    plt.show()

def changeTaskId(tid):
    fname = 'input.br.in'
    for line in fileinput.input(fname, inplace=True):
        words = line.split()
        if len(words) == 3 and words[0] == "taskID":
            line = line.replace(words[2], str(tid))
        print line,
    return
            

def run():
    startTime = time.time()
    taskId = range(100)
    for tid in taskId:
        changeTaskId(tid)
        os.system('make run')
    endTime = time.time()
    print "Total Running Time: ", endTime - startTime

def main():
    # run()
    plotFig()

if __name__ == "__main__":
    main()

