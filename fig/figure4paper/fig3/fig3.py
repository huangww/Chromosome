import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure(0,figsize=(8,3))
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 12 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text',usetex = True)

def raw_mean(i, N, T):
    mu = (N+1)/2.0
    raw_mean = T*np.log(((i-mu)/T)/(np.sinh((i-mu)/T)))
    return raw_mean

def z_mean(i, N, T):
    z_mean = raw_mean(i,N,T)-raw_mean(0,N,T)
    return z_mean

def integral_xy_var(i, N, T):
    mu = (N+1)/2.0
    mean_z = np.cosh((mu-i)/T)/np.sinh((mu-i)/T)-T/(mu-i)
    var_z = -(1.0/np.sinh((mu-i)/T))**2+(T/(mu-i))**2
    second_moment_z = var_z + mean_z*mean_z
    integral_xy_var = 1 - second_moment_z
    return integral_xy_var

def xy_raw_var(m, n, N, T):
    x = np.linspace(0,N,N)
    integral = integral_xy_var(x,N,T)
    xy_raw_var = sum(integral[m:n])
    return xy_raw_var

def z_raw_var(i, N, T):
    mu = (N+1)/2.0
    z_raw_var = T*(np.cosh((i-mu)/T)/np.sinh((i-mu)/T)-1.0/((i-mu)/T))
    return z_raw_var

def xy_var(i, N, T):
    xy_var = xy_raw_var(0,i,N,T)*xy_raw_var(i,N,N,T)/xy_raw_var(0,N,N,T)
    return xy_var

def z_var(i, N, T):
    z_var = (z_raw_var(i,N,T)-z_raw_var(0,N,T))*(z_raw_var(N,N,T)-z_raw_var(i,N,T))/(z_raw_var(N,N,T)-z_raw_var(0,N,T))
    return z_var

N = 50
index = np.linspace(0,N,N)
sp1 = plt.subplot(121)
sp2 = plt.subplot(122)
for T in [1,5,10,50]:
    data = np.loadtxt('MD_N50_T'+str(T)+'_pseudo_xyz_varxyz.dat')
    line, = sp1.plot(index,z_mean(index,N,T))
    sp1.plot(index[::2],data[0,::2],line.get_color()+'o')
    sp1.fill_between(index,
            z_mean(index,N,T)-np.sqrt(z_var(index,N,T)),
            z_mean(index,N,T)+np.sqrt(z_var(index,N,T)),
            facecolor = 'green',
            alpha = 0.5)
# sp1.set_xlabel('Bead index')
# sp1.set_ylabel('Average $z$')

N = 500
index = np.linspace(0,N,N)
T = np.linspace(0.1,1000,10000)
var_xy = np.zeros(len(T))
# for l in [50,100,200,250]:
for l in [250]:
    i = 0
    for t in T:
        var_xy[i] = xy_var(l,N,t)
        i = i + 1
    var_z = z_var(l,N,T)
    sp2.semilogx(T,np.sqrt(2.0*(var_z+var_xy)))
# sp2.set_yticks([0, 1, 2, 3, 4])
# sp2.set_xlim([0,1000])
# sp2.set_xlabel('Effective Temperature')
# sp2.set_ylabel('Distance between loci pair')
fig.show()
fig.savefig('fig3.pdf')
