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

N = 500
T = '200'
dataDir = 'data/'
index = np.linspace(0,N,N)
sp1 = plt.subplot(121)
sp2 = plt.subplot(122)
# for T in [1,5,10,50]:
dataFile = dataDir + 'MD_N'+str(N)+'_T'+str(T)+'_xyz_varxyz.dat'
data = np.loadtxt(dataFile)
line, = sp1.plot(index,z_mean(index,N,int(T)))
sp1.plot(index[::5],data[::5,0],line.get_color()+'o')
# var_z = np.zeros(len(index))
# for i in range(N):
    # var_z[i] = z_var(i,N,int(T)
# line, = sp2.plot(index,var_z)
line, = sp2.plot(index,z_var(index,N,int(T)))
sp2.plot(index[::5],data[::5,3],line.get_color()+'o')

plt.show()
