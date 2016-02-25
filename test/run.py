import numpy as np
import fileinput
import os 
import itertools as it
import time

Teff = np.linspace(10, 100, 91)
Nbead = [100, 500, 1000]
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
            line = line.replace(words[2], str(T/2))
        if len(words) == 3 and words[0] == "tEnd":
            line = line.replace(words[2], str(T*1e6))
        print line,

if __name__ == "__main__":
    startTime = time.time()
    for (N,T) in it.product(Nbead, Teff):
        # print N, T
        ChangeInput(N, T)
        os.system('make run')
    endTime = time.time()
    print "Total Running Time: ", endTime - startTime
