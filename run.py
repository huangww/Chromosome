#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import fileinput
import os 
import itertools as it
import time
import glob
import script.acfFit as af

# Teff = np.linspace(0.1, 9.9, 99)
# Teff = np.concatenate([np.linspace(1, 9.5, 18), np.linspace(10, 95, 18), np.linspace(100, 1000, 19)])
Teff = [100, 200] 
Nbead = [100]
# Nbead = [10, 20, 30, 50, 70, 100, 200, 300, 500]
# Teff = [0] 
fname = 'input.in'

def ChangeInput(N, T):
    for line in fileinput.input(fname, inplace=True):
    # for line in open(fname):
        words = line.split()
        if len(words) == 3 and words[0] == "nSite":
            line = line.replace(words[2], str(N))
        if len(words) == 3 and words[0] == "nPar":
            line = line.replace(words[2], str(N/2))
        if len(words) == 3 and words[0] == "tempEff":
            line = line.replace(words[2], str(T))
        if len(words) == 3 and words[0] == "dt":
            dt = min(T/10., 100)
            # dt = min(N/100., 100)
            # line = line.replace(words[2], str(T/2))
            line = line.replace(words[2], str(dt))
        if len(words) == 3 and words[0] == "tEnd":
            # line = line.replace(words[2], str(T*1e6))
            line = line.replace(words[2], str(dt*1e6))
        print line,

def GetDataTau(N, T):
    fname = 'data/Asep_rg_N%g_T%g_*.dat'%(N, T)
    fname = glob.glob(fname)[0]
    print fname
    data = np.loadtxt(fname)
    t = data[:, 0]
    tau = [T, af.GetTau(t, data[:,1]), af.GetTau(t, data[:,-1])]
    return tau

def GetTaskId(fname):
    fname = glob.glob(fname)[0]
    split1 = fname.split('_')
    split2 = split1[-1].split('.')
    return split2[0]
    

def main():
    startTime = time.time()
    tau = []
    N = Nbead[0]
    for T in Teff:
    # for (N,T) in it.product(Nbead, Teff):
        # print N, T
        ChangeInput(N, T)
        os.system('make run')
        tau.append(GetDataTau(N, T))
    fname = 'data/Asep_rg_N%g_T%g_*.dat'%(N, T)
    fname = 'data/Asep_tau_N%g_'%N + GetTaskId(fname) + '.dat'
    np.savetxt(fname, tau)
    endTime = time.time()
    print "Total Running Time: ", endTime - startTime

if __name__ == "__main__":
    main()
