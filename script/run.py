#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import fileinput
import os 
import sys
import subprocess
import itertools as it
import time
import script.acfFit as af
import glob

# Teff = np.concatenate([np.linspace(10, 90, 9), np.linspace(100, 1000, 10)])
# Teff = np.concatenate([np.linspace(1, 9.5, 18), np.linspace(10, 49, 40), np.linspace(50,190,15), np.linspace(200, 1000, 17)])
Teff = range(1, 30)
Nbead = [100]
infile = 'input.br.in'

def ChangeInput(N, T):
    if T < 25:
        dt = 1e-4
    elif T >= 25 and T < 80:
        dt = 2e-5
    else:
        dt = 1e-5
    outStep = 1e4
    fExternal = 1
    for line in fileinput.input(infile, inplace=True):
        words = line.split()
        if len(words) == 3 and words[0] == "nBead":
            line = line.replace(words[2], str(N))
        if len(words) == 3 and words[0] == "fExternal":
            line = line.replace(words[2], str(fExternal))
        if len(words) == 3 and words[0] == "temperature":
            line = line.replace(words[2], str(T))
        if len(words) == 3 and words[0] == "dt":
            line = line.replace(words[2], str(dt))
        if len(words) == 3 and words[0] == "outputStep":
            line = line.replace(words[2], str(outStep))
        if len(words) == 3 and words[0] == "tEnd":
            line = line.replace(words[2], str(dt*1e8))
        print line,

def Changedt(dt):
    for line in fileinput.input(infile, inplace=True):
        words = line.split()
        if len(words) == 3 and words[0] == "dt":
            line = line.replace(words[2], str(dt))
        if len(words) == 3 and words[0] == "outputStep":
            line = line.replace(words[2], str(1e5))
        if len(words) == 3 and words[0] == "tEnd":
            line = line.replace(words[2], str(dt*1e4))
        print line,
    return

def Finddt():
    dtUp = 1e-3
    dtLow = 1e-6
    dt = dtUp
    nullFile = open(os.devnull, 'w')
    while dtUp - dtLow > 1e-6:
        Changedt(dt)
        try:
            # exitStatus = subprocess.call(["make","run"])
            exitStatus = subprocess.call(["make","run"],stdout=nullFile, stderr=subprocess.STDOUT)
            if exitStatus == 0:
                dtLow = dt
            else:
                dtUp = dt
        except Exception:
            pass
        dt = (dtUp + dtLow) / 2.
        # print dtLow, dt, dtUp
    dtString = ("%.1e" % dt).split('.')
    dt = dtString[0] + dtString[1][1:]
    return dt
        

def TreatData(N, T):
    fname = 'data/BeadRod_rg_N%g_F%g_T%g_*.dat'%(N, 1, T)
    fname = glob.glob(fname)[0] 
    print fname
    data = np.loadtxt(fname)
    t = data[:,0]
    x = data[:,1]
    tau = [T, af.GetTau(t, x)]
    # fname = 'data/r_N%g_T%g_*.dat'%(N, T)
    # fname = glob.glob(fname)[0] 
    # print fname
    # data = np.loadtxt(fname).reshape([-1,N,3])
    # tagBead = [1, 10, 20, 30, 40, 50]
    # for i in tagBead:
    #     x = data[:, i, 0]
    #     y = data[:, i, 1]
    #     z = data[:, i, 2]
    #     tau.append(af.GetTau(t, x))
    #     tau.append(af.GetTau(t, y))
    #     tau.append(af.GetTau(t, z))
    fname = fname.replace('rg', 'tau')
    np.savetxt(fname, tau)
    
if __name__ == "__main__":
    startTime = time.time()
    for (N,T) in it.product(Nbead, Teff):
        ChangeInput(N, T)
        print T, Finddt()
        # try:
        #     exitStatus = subprocess.call(["make","run"])
        # except Exception:
        #     pass
        # if exitStatus == 0:
        #     TreatData(N, T)
    endTime = time.time()
    print "Total Running Time: ", endTime - startTime
