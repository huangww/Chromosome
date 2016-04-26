#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import fileinput
import os 
import itertools as it
import time

# Teff = np.linspace(0.1, 9.9, 99)
Teff = np.concatenate([np.linspace(1, 9.5, 18), np.linspace(10, 95, 18), np.linspace(100, 1000, 19)])
# Teff = [2000] 
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

if __name__ == "__main__":
    startTime = time.time()
    for (N,T) in it.product(Nbead, Teff):
        # print N, T
        ChangeInput(N, T)
        os.system('make run')
    endTime = time.time()
    print "Total Running Time: ", endTime - startTime
