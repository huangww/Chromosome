#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import fileinput
import os 
import itertools as it
import time
import glob
import script.acfFit as af
import statsmodels.tsa.stattools as ss

fname = 'input.in'

def ChangeInput(T, tid):
    for line in fileinput.input(fname, inplace=True):
        words = line.split()
        if len(words) == 3 and words[0] == "tempEff":
            line = line.replace(words[2], str(T))
        if len(words) == 3 and words[0] == "taskID":
            line = line.replace(words[2], str(tid))
        if len(words) == 3 and words[0] == "dt":
            line = line.replace(words[2], str(T))
        if len(words) == 3 and words[0] == "tEnd":
            line = line.replace(words[2], str(T*1e7))
        print line,

def GetACF(fname):
    data = np.loadtxt(fname)
    t = data[:,0]
    x = data[:,-1]
    acf= ss.acf(x, nlags=len(x), fft=True)
    return t[:len(acf)], acf

def main():
    startTime = time.time()
    N = 100
    Teff = [20] 
    taskID = [5]
    for (T, tid) in it.product(Teff, taskID):
        ChangeInput(T, tid)
        os.system('make run')
        fname = 'data/Asep_rg_N100_T%g_%g.dat'%(T, tid)
        t, acf = GetACF(fname)
        fname = fname.replace('rg', 'acf')
        print fname
        np.savetxt(fname, [t, acf])
    endTime = time.time()
    print "Total Running Time: ", endTime - startTime

if __name__ == "__main__":
    main()
