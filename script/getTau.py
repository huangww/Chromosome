"""
This script treat the trajectory data and get information like relaxation time scale
"""
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import acfFit as af
import itertools as it



N = 100
# Teff = range(1, 10) + range(10, 20, 2) + range(20, 30, 5) + range(30, 40, 10) + \
#         range(40, 100, 20) + range(100, 200, 50) + range(200, 400, 100) + \
#         range(400, 1000, 200) + [1000]
# Teff = [1, 5]
Teff = range(2, 32) + range(35, 155, 5)
Teff.remove(20)
taskId = range(1, 101)

def GetT(fname):
    paras = fname.split('_')
    T = [float(s[1:]) for s in paras if s[0]=='T']
    return T[0]

def ExtractTau(data):
    t = data[:,0]
    tau = []
    for i in range(1, len(data[0])):
        x = data[:,i]
        tau.append(af.GetTau(t, x))
    return np.array(tau)


def MergeTau(N, T):
    fname = dataDir + 'tau_N%g_T%g_*.dat'%(N, T)
    fileList = glob.glob(fname)
    tauList = []
    for f in fileList:
        tau = np.loadtxt(f)
        tau[tau>1e4] = np.nan
        tau[tau<0] = np.nan
        tauList.append(tau)
    tau = np.nanmean(tauList, axis=0)
    tauVar = np.nanvar(tauList, axis=0)
    return tau, tauVar

def PlotFig(tauData):
    ax = plt.figure().add_subplot(111)
    for i in range(1, len(tauData[0])):
        ax.plot(tauData[:,0], tauData[:,i], '-o')
    labels = [r'$\tau_{R_g}$',  r'$\tau_{1}$',
            r'$\tau_{10}$', r'$\tau_{20}$',
            r'$\tau_{30}$', r'$\tau_{40}$',
            r'$\tau_{50}$']
    ax.legend(labels, ncol=4)
    ax.set_xlabel(r'$1/F$')
    ax.set_ylabel(r'$<\tau>$')
    ax.set_title(r'$N = 100$')
    plt.savefig('fig/tau_BS_N100.pdf')

def PlotFig3D(tauData):
    ax = plt.figure().add_subplot(111)
    for i in range(1, len(tauData[0])):
        ax.plot(tauData[:,0], tauData[:,i], '-o')
    labels = [r'$\tau_{R_g}$',  r'$\tau_{1}$',
            r'$\tau_{10}$', r'$\tau_{20}$',
            r'$\tau_{30}$', r'$\tau_{40}$',
            r'$\tau_{50}$']
    ax.legend(labels, ncol=4)
    ax.set_xlabel(r'$1/F$')
    ax.set_ylabel(r'$<\tau>$')
    ax.set_title(r'$N = 100$')
    plt.savefig('fig/tau_BS_N100.pdf')


def GetTauFromTauData():
    tauData = []
    tauVarData = []
    for T in Teff:
        tau,tauVar = MergeTau(N, T)
        tau = np.insert(tau, 0, T)
        tauVar = np.insert(tauVar, 0, T)
        # print tau
        tauData.append(tau)  
        tauVarData.append(tauVar)  
    np.savetxt('data/tau1D_N100.dat', tauData)
    np.savetxt('data/tauVar1D_N100.dat', tauVarData)
    tauData = np.loadtxt('data/tau1D_N100.dat')
    PlotFig(tauData)

def GetTauFromRawData():
    tauData = []
    for T in Teff:
        fname = 'beadRodN100/BeadRod_rd_N100_F0.1_T%g_*.dat'%T
        fList = [f for f in glob.iglob(fname) if os.path.getsize(f) > 1e5]
        if fList == []:
            pass
        else:
            tauList = []
            for f in fList:
                print f
                data = np.loadtxt(f)
                tau = ExtractTau(data)
                tau[tau>1e3] = np.nan
                tau[tau<0] = np.nan
                tauList.append(tau)
            tau = np.nanmean(tauList, axis=0)
            tau = np.insert(tau, 0, T)
            tauData.append(tau.tolist())
    fname = 'data/tauNoise_N100_F0.1.dat'
    np.savetxt(fname, np.array(tauData))
    tauData = np.loadtxt(fname)
    PlotFig(tauData)

def GetTau(fname):
    fnameList = glob.glob(fname)
    tau = []
    for f in fnameList:
        print f
        t = np.loadtxt(f.replace('rd', 'rg'))[:,0]
        x = np.loadtxt(f)[:,0]
        tau0 = af.GetTau(t, x)
        # if tau0 > 200 and tau0 < 500:
        if not(np.isnan(tau0)) and tau0<600:
            tau.append(tau0)
    return np.mean(np.array(tau))


def GetTauMean():
    tau = []
    for T in Teff:
        # fname = 'BR-N100/rd_N100_T%g_*.dat'%T
        fname = 'beadRodN100Teff/BeadRod_rd_N100_T%g_*.dat'%T
        tau.append([T,GetTau(fname)])
    print tau
    np.savetxt('BeadRod_tau_N100_Teff.dat', tau)


def main():
    GetTauFromRawData()
    # GetTauFromTauData()
    # GetTauMean()

if __name__ == "__main__":
    main()        


