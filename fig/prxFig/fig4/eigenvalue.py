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


def Demo():            
    fig = plt.figure(0,figsize=(5, 3.6))
    plt.rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'serif'  : 'Helvetica',
            'weight' : 'normal',
            'size'   : 15 }
    plt.rc('lines', lw=2)
    plt.rc('font', **font)
    # fig, (ax1, ax2) = plt.subplots(1,2, gridspec_kw = 
    #         {'width_ratios':[1, 3]})
    fig.subplots_adjust(left=0.13, right =0.95,\
            bottom=0.13, top =0.95, wspace=0.15)

    ax1 = plt.subplot2grid((1, 3), (0, 0))
    a, b, L = (1, 2, 10)
    M = TransitionMartrix(a, b, 1, L)
    ev1 = eigvals(M)/(a+b)
    for x in ev1:
        ax1.scatter(1, x, marker='_', color='r')
    M = TransitionMartrix(a, b, 2, L)
    ev2 = eigvals(M)/(a+b)
    for x in ev2:
        if np.isclose(ev1, x, rtol=1e-3).any():
            ax1.scatter(2, x, marker='_', color='r')
        else:
            ax1.scatter(2, x, marker='_', color='b')
    M = TransitionMartrix(a, b, 3, L)
    ev3 = eigvals(M)/(a+b)
    for x in ev3:
        if np.isclose(ev1, x, rtol=1e-3).any():
            ax1.scatter(3, x, marker='_', color='r')
        elif np.isclose(ev2, x, rtol=1e-3).any():
            ax1.scatter(3, x, marker='_', color='b')
        else:
            ax1.scatter(3, x, marker='_', color='g')
    ax1.set_xlabel('Number of particles')
    ax1.set_xticks([1, 2, 3])
    ax1.set_ylabel('Eigenvlues')
    ax1.set_ylim(ymax=1)

    ax2 = plt.subplot2grid((1,3), (0, 1), colspan=2)
    ev2 = np.sort(ev2)[::-1]
    for i,x in zip(range(len(ev2)), ev2):
        ax2.scatter(i, x, s=30, marker='s', color='b')
        ax2.scatter(i, x, s=20, marker='*', color='r')
    ax2.set_xlabel(r'$k$')
    # ax2.set_ylabel('Eigenvalues')
    ax2.legend(['Matrix diagonalization', 'Bethe-ansatz'], loc='lower left', fontsize=14)
    ax2.set_ylim(ax1.get_ylim())
    ax2.set_xlim([0, 50])
    ax2.set_yticklabels([])
    
    plt.show()

def main():
    # EigenValues()
    # ExtraEigenVector()
    # MinimalExample()
    Demo()



if __name__ == "__main__":
    main()
