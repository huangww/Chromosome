import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import comb
import itertools as it
from numpy.linalg import *

def TransitionMartrixPeriodic(a, b, N, L):
    nDim = int(comb(L, N))
    M = np.zeros([nDim, nDim])
    pos = range(L)
    for i,x in zip(range(nDim),it.combinations(pos, N)):
        for j,y in zip(range(nDim),it.combinations(pos, N)):
            diffxy = np.array(x) - np.array(y)
            if sum(abs(diffxy)) == 1:
                if sum(diffxy) == 1:
                    M[i,j] = a
                if sum(diffxy) == -1:
                    M[i,j] = b
            if (y[-1]==(L-1) and x[0]==0 and x[1:]==y[:-1]):
                M[i,j] = a
            if (y[0]==0 and x[-1]==(L-1) and x[:-1]==y[1:]):
                M[i,j] = b
    for k in range(len(M)):
        M[k,k] = -sum(M[:,k])
    return M


def TransitionMartrix(a, b, N, L):
    nDim = int(comb(L, N))
    M = np.zeros([nDim, nDim])
    pos = range(L)
    for i,x in zip(range(nDim),it.combinations(pos, N)):
        for j,y in zip(range(nDim),it.combinations(pos, N)):
            diffxy = np.array(x) - np.array(y)
            if sum(abs(diffxy)) == 1:
                if sum(diffxy) == 1:
                    M[i,j] = a
                if sum(diffxy) == -1:
                    M[i,j] = b
    for k in range(len(M)):
        M[k,k] = -sum(M[:,k])
    return M

def EigenValues():            
    L = 10 
    a = 1
    b = 2
    ax = plt.figure(0).add_subplot(111)
    ax.axhline(y=-(a+b)+2*np.sqrt(a*b)*np.cos(np.pi/L), color='r')
    # M = TransitionMartrix(a, b, 2, L)
    # print M
    for i in range(1, L):
        M = TransitionMartrix(a, b, i, L)
        ev = eigvals(M)
        for x in ev:
            ax.scatter(i, x,marker='_')
    ax.set_xlabel('Number of particles')
    ax.set_ylabel('Eigenvlues')
    # ax.set_ylim([-5,2])
    plt.show()

def ExtraEigenVector():
    L = 10
    a = 1
    b = 2
    pos = range(L)
    M = TransitionMartrix(a, b, 1, L)
    vals1, vectors1 = eig(M)
    N = 2
    M = TransitionMartrix(a, b, N, L)
    vals2, vectors2 = eig(M)
    for x2 in vals2:
        if  (np.any(abs(x2-vals1)<1e-6)):
            print x2
            v = vectors2[np.where(vals2==x2)]
            if abs(x2) > 1e-4:
                plt.plot(v[0,:L],'-o')
                break
            # count = 0
            # data = []
            # for config in it.combinations(pos, N):
            #     a =  [i for i in config]
            #     a.append(v[0,count])
            #     data.append(a)
            #     count = count + 1
            # data = np.array(data)
            # plt.plot(data[:,0], data[:,1], data[:,-1])
    plt.show()

def MinimalExample():
    N, L = (2, 4)
    a, b = (2, 1)
    M = TransitionMartrixPeriodic(a, b, N, L)
    print M
    ev = eigvals(M)
    print ev


def main():
    # EigenValues()
    # ExtraEigenVector()
    MinimalExample()


if __name__ == "__main__":
    main()
