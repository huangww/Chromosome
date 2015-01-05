import numpy as np

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

def r_var(i, N, T):
    r_var = xy_var(i,N,int(T)) + z_var(i,N,int(T))
    return r_var

def sum_r_var(N, T):
    var_r = np.zeros(N)
    for i in range(N):
        var_r[i] = r_var(i,N,int(T))
    sum_r_var = sum(var_r)
    return sum_r_var

def Gyration(N, T):
    Rgs = sum_r_var(N, T)/(N*N)
    return Rgs
    
