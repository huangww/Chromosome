#!/usr/bin/env python
# encoding: utf-8
import numpy as np

def mean_n_1D(i, N, T):
    Teff = float(T)/2
    mu = (N+1)/2.0
    mean_n_1D = 1.0/(1+np.exp((i-mu)/Teff))
    return mean_n_1D

def mean_zi_1D(i, N, T):
    Teff = float(T)
    mean_zi_1D =Teff*np.log(((1+np.exp(N/Teff))/
            (np.exp(i/Teff)+np.exp((N-i)/Teff))))
    return mean_zi_1D

def mean_z_1D_sum(N, T):
    Teff = float(T)
    mean_z_1D_sum = np.zeros(N)
    for i in range(N):
        mean_z_1D_sum[i] = 0
        for i0 in range(i):
            mean_z_1D_sum[i] += (2*mean_n_1D(i0+1, N, Teff)-1)
    return mean_z_1D_sum

def mean_z_1D(N, T):
    Teff = float(T)
    index = np.arange(N)
    mean_z_1D = Teff*np.log(((1+np.exp(N/Teff))/
            (np.exp(index/Teff)+np.exp((N-index)/Teff))))
    return mean_z_1D

def var_z_1D(N, T):
    Teff = float(T)
    index = np.arange(N)
    var_z_1D = Teff*np.sinh((N-index)/Teff) \
            *np.sinh(index/Teff)/(np.sinh(N/Teff)* \
            np.cosh((N-2*index)/(2*Teff))**2)
    return var_z_1D

def mean_e_z(i, N, T):
    mu = (N+1)/2.0
    x = (mu - i)/float(T)
    if x == 0:
        mean_e_z = 0
    else:
        mean_e_z = 1.0/np.tanh(x) - 1.0/x
    return mean_e_z

def var_e_z(i, N, T):
    mu = (N+1)/2.0
    x = (mu - i)/float(T)
    var_e_z = 1/x**2 - 1/np.sinh(x)**2 
    return var_e_z

def mean_zi(i, N, T):
    mean_z = 0
    for i0 in range(i):
        mean_z += mean_e_z(i0+1, N, T)
    return mean_z

def mean_z(N, T):
    z_mean = np.zeros(N)
    for i in range(N):
        z_mean[i] = 0
        for i0 in range(i):
            z_mean[i] += mean_e_z(i0+1, N, T)
    return z_mean

def meanz(N, T):
    zmean = np.zeros(N)
    Teff = float(T)
    for i in range(N):
        mu = i - (N+1)/2.0
        zmean[i] = Teff * np.log((2*mu*np.sinh(N/(2*Teff)))/
            (N*np.sinh(mu/Teff)))
    return zmean

def integral_xy_var(i, N, T):
    mu = (N+1)/2.0
    mean_e_z = 1.0/np.tanh((mu-i)/T)-T/(mu-i)
    var_e_z = -(1.0/np.sinh((mu-i)/T))**2+(T/(mu-i))**2
    second_moment_e_z = var_e_z + mean_e_z*mean_e_z
    integral_xy_var = 1 - second_moment_e_z
    return integral_xy_var

def xy_var(m, n, N, T):
    x = np.linspace(0,N,N)
    integral = integral_xy_var(x,N,T)
    xy_var = sum(integral[m:n])
    return xy_var

# need to check again
def z_var(m, n, N, T):
    # mu = (N+1)/2.0
    # z_var = T*((1.0/np.tanh((n-mu)/T)-1.0/((n-mu)/T))-(1.0/np.tanh((m-mu)/T)-1.0/((m-mu)/T)))
    z_var = 0
    for i in range(m,n):
        z_var += var_e_z(i, N, T)
    return z_var

def var_z(N, T):
    var_z = np.zeros(N+1)
    for i in range(N+1):
        v0s = z_var(0,i,N,T)
        vst = z_var(i,N,N,T)
        v0t = v0s + vst
        # var_z[i] = 2.0*v0s*vst/v0t
        var_z[i] = v0s*vst/v0t
    return var_z

def var_xy(N, T):
    var_xy = np.zeros(N+1)
    for i in range(N+1):
        v0s = xy_var(0,i,N,T)
        vst = xy_var(i,N,N,T)
        v0t = v0s + vst
        # var_xy[i] = 2.0*v0s*vst/v0t
        var_xy[i] = v0s*vst/v0t
    return var_xy

def var_xy_cm(cm, N, T):
    var_xy_cm = np.zeros(N+1)
    for i in range(N+1):
        if i<cm:
            v0s = xy_var(0,i,N,T)
            vst = xy_var(i,cm,N,T)
            v0t = v0s + vst
            var_xy_cm[i] = v0s*vst/v0t
        else:
            v0s = xy_var(cm,i,N,T)
            vst = xy_var(i,N,N,T)
            v0t = v0s + vst
            var_xy_cm[i] = v0s*vst/v0t
    return var_xy_cm

def var_z_cm(cm, N, T):
    var_z_cm = np.zeros(N+1)
    for i in range(N+1):
        if i<cm:
            v0s = z_var(0,i,N,T)
            vst = z_var(i,cm,N,T)
            v0t = v0s + vst
            var_z_cm[i] = 2*v0s*vst/v0t
        else:
            v0s = z_var(cm,i,N,T)
            vst = z_var(i,N,N,T)
            v0t = v0s + vst
            var_z_cm[i] = 2*v0s*vst/v0t
    return var_z_cm

