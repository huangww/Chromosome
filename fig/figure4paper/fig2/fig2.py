import numpy as np 
import matplotlib.pyplot as plt

def x_mean(x, N, T):
    x_mean = T*np.log((1.0+np.exp(N/float(T)))/(np.exp(x/float(T))+np.exp((N-x)/float(T))))
    return x_mean

def raw_var(x, N, T):
    raw_var = 2.0*T/(1.0+np.exp(2.0*(x-N/2)/float(T)))
    return raw_var

def x_var(x, N, T):
    x_var = (raw_var(0,N,T)-raw_var(x,N,T))*(raw_var(x,N,T)-raw_var(N,N,T))/(raw_var(0,N,T)-raw_var(N,N,T))
    return x_var
    
fig = plt.figure(0,figsize=(8,3))
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 12 }
plt.rc('lines', lw=2)
plt.rc('font', **font)

N = 500
x = np.linspace(1,N,N)
sp1 = plt.subplot(121)
sp2 = plt.subplot(122)
# for T in [100,300,1000,3000]:
for T in [1,100,200,5000]:
    # data = np.loadtxt('visualization_T'+str(T)+'_N500.txt')
    # meanx = data.mean(axis=0)
    # varx = data.var(axis=0)
    line, = sp1.plot(x,x_mean(x,N,T))
    # sp1.plot(x[::50],meanx[::5],line.get_color()+'o')
    sp2.plot(x,x_var(x,N,T),line.get_color())
    # sp2.plot(x[::50],varx[::5],line.get_color()+'o')
#sp1.set_xlabel('Bead index')
#sp1.set_ylabel('Average position')
#sp2.set_xlabel('Bead index')
#sp2.set_ylabel('Varriance')
# sp1.set_yticks([100, 200, 300, 400, 500])
# fig.show()
plt.show()
# fig.savefig('fig2_bc.eps')

