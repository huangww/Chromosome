import numpy as np
import matplotlib.pyplot as plt


fig = plt.figure(0,figsize=(4,3))
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 12 }
plt.rc('lines', lw=2)
plt.rc('font', **font)

def integral_xy_var(i, N, T):
    mu = (N+1)/2.0
    mean_z = np.cosh((mu-i)/T)/np.sinh((mu-i)/T)-T/(mu-i)
    var_z = -(1.0/np.sinh((mu-i)/T))**2+(T/(mu-i))**2
    second_moment_z = var_z + mean_z*mean_z
    integral_xy_var = 1 - second_moment_z
    return integral_xy_var

def xy_var(m, n, N, T):
    x = np.linspace(0,N,N)
    integral = integral_xy_var(x,N,T)
    xy_var = sum(integral[m:n])
    return xy_var

def z_var(m, n, N, T):
    mu = (N+1)/2.0
    z_var = T*((np.cosh((n-mu)/T)/np.sinh((n-mu)/T)-1.0/((n-mu)/T))-(np.cosh((m-mu)/T)/np.sinh((m-mu)/T)-1.0/((m-mu)/T)))
    return z_var

N = 50
centromere_location = 19
x = np.linspace(1,N,N)
var_z = np.zeros(N)
var_xy = np.zeros(N)
for T in [1,10,50]:
    for i in range(N):
        if i<centromere_location:
            v0s = z_var(0,i,N,T)
            vst = z_var(i,centromere_location,N,T)
            v0t = v0s + vst
            var_z[i] = 2*v0s*vst/v0t
            #v0s = xy_var(0,i,N,T)
            #vst = xy_var(i,centromere_location,N,T)
            #v0t = v0s + vst
            #var_xy[i] = v0s*vst/v0t
        else:
            v0s = z_var(centromere_location,i,N,T)
            vst = z_var(i,N,N,T)
            v0t = v0s + vst
            var_z[i] = 2*v0s*vst/v0t
            #v0s = xy_var(centromere_location,i,N,T)
            #vst = xy_var(i,N,N,T)
            #v0t = v0s + vst
            #var_xy[i] = v0s*vst/v0t
    line, = plt.plot(x,var_z)
    filename = 'MD_N98_T'+str(T)+'_pseudo_xyz_varxyz.dat'
    data = np.loadtxt(filename)
    plt.plot(x[::2],data[3,::2],line.get_color()+'o')
T = 5000
for i in range(N):
    v0s = z_var(0,i,N,T)
    vst = z_var(i,N,N,T)
    v0t = v0s + vst
    var_z[i] = 2.0*v0s*vst/v0t
plt.plot(x,var_z,'k--')
plt.ylim([0,10])

fig.savefig('fig4.pdf')
fig.show()
