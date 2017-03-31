#!/usr/bin/env python
# encoding: utf-8
import numpy as np

# Theory for the 1D pinned polymer loop
def mean_n_1D(i, N, T):
    Teff = float(T)/2
    mu = (N-1)/2.0
    mean_n_1D = 1.0/(1+np.exp((i-mu)/Teff))
    return mean_n_1D

def var_n_1D(i, N, T):
    p = mean_n_1D(i, N, T)
    var_n_1D = p * (1 - p)
    return var_n_1D

def mean_zi_1D(i, N, T):
    Teff = float(T)
    mean_zi_1D = 0
    for i0 in range(i):
        mean_zi_1D += (2*mean_n_1D(i0, N, Teff)-1)
    return mean_zi_1D

def z_var_1D(m, n, N, T):
    z_var_1D = 0
    for i in range(m,n):
        z_var_1D += 4*var_n_1D(i, N, T)
    return z_var_1D

def var_zi_1D(i, N, T):
    Teff = float(T)
    var_zi_1D = 0
    for i0 in range(i):
        v0s = z_var_1D(0,i0,N,T)
        vst = z_var_1D(i0,N,N,T)
        v0t = v0s + vst
        var_zi_1D = v0s*vst/v0t
    return var_zi_1D

def mean_z_1D_sum(N, T):
    Teff = float(T)
    mean_z_1D_sum = np.zeros(N)
    for i in range(N):
        mean_z_1D_sum[i] = 0
        for i0 in range(i):
            mean_z_1D_sum[i] += (2*mean_n_1D(i0, N, Teff)-1)
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

# Theory for the 3D pinned polymer loop
def mean_e_z(i, N, T):
    mu = (N-1)/2.0
    x = (mu - i)/float(T)
    mean_e_z = 1.0/np.tanh(x) - 1.0/x
    # mean_e_z = 0
    return mean_e_z

def var_e_z(i, N, T):
    mu = (N-1)/2.0
    x = (mu - i)/float(T)
    var_e_z = 1/x**2 - 1/np.sinh(x)**2 
    return var_e_z

def mean_zi(i, N, T):
    mean_z = 0
    for i0 in range(i):
        mean_z += mean_e_z(i0, N, T)
    return mean_z

def mean_z(N, T):
    mean_z = np.zeros(N)
    for i in range(N):
        mean_z[i] = 0
        for i0 in range(i):
            mean_z[i] += mean_e_z(i0, N, T)
    return mean_z

def meanz(N, T):
    zmean = np.zeros(N)
    Teff = float(T)
    for i in range(N):
        mu = i - (N-1)/2.0
        zmean[i] = Teff * np.log((2*mu*np.sinh(N/(2*Teff)))/
            (N*np.sinh(mu/Teff)))
    return zmean

def integral_xy_var(i, N, T):
    mu = (N-1)/2.0
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
    # mu = (N-1)/2.0
    # z_var = T*((1.0/np.tanh((n-mu)/T)-1.0/((n-mu)/T))-(1.0/np.tanh((m-mu)/T)-1.0/((m-mu)/T)))
    z_var = 0
    for i in range(m,n):
        z_var += var_e_z(i, N, T)
    return z_var

def var_zi(i, N, T):
    v0s = z_var(0,i,N,T)
    vst = z_var(i,N,N,T)
    v0t = v0s + vst
    var_zi = v0s*vst/v0t
    return var_zi

def var_z(N, T):
    var_z = np.zeros(N)
    for i in range(N):
        v0s = z_var(0,i,N,T)
        vst = z_var(i,N,N,T)
        v0t = v0s + vst
        var_z[i] = v0s*vst/v0t
    return var_z

def var_xyi(i, N, T):
    v0s = xy_var(0,i,N,T)
    vst = xy_var(i,N,N,T)
    v0t = v0s + vst
    # var_xy[i] = 2.0*v0s*vst/v0t
    var_xyi = v0s*vst/v0t
    return var_xyi

def var_xy(N, T):
    var_xy = np.zeros(N)
    for i in range(N):
        v0s = xy_var(0,i,N,T)
        vst = xy_var(i,N,N,T)
        v0t = v0s + vst
        # var_xy[i] = 2.0*v0s*vst/v0t
        var_xy[i] = v0s*vst/v0t
    return var_xy

def var_ri(i, N, T):
    return var_xyi(i, N, T) + var_zi(i, N, T)

def var_r(N, T):
    return var_xy(N, T) + var_z(N, T)

def var_xy_cm(cm, N, T):
    var_xy_cm = np.zeros(N)
    for i in range(N):
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
    var_z_cm = np.zeros(N)
    for i in range(N):
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

# Particle density function
def particle_density_z(z, N, T):
    zmean = mean_z(N, T)
    zvar = var_z(N, T)
    pz = np.zeros([len(zmean),len(z)])
    import matplotlib.mlab as mlab
    for i in range(len(zmean)):
        pz[i] = mlab.normpdf(z, zmean[i], np.sqrt(zvar[i]))
    pdf = np.nanmean(pz, axis=0)
    return pdf
    

def particle_density_x(x, N, T):
    xvar = var_xy(N, T)/2.
    px = np.zeros([len(xvar),len(x)])
    import matplotlib.mlab as mlab
    for i in range(len(xvar)):
        px[i] = mlab.normpdf(x, 0, np.sqrt(xvar[i]))
    pdf = np.nanmean(px, axis=0)
    return pdf

# Theory for the Rouse model
def rouse_mean(N, T):
    f = 1./float(T)
    k_H = 3
    z_mean = np.zeros(N)
    l = np.arange(1, (N+1)/2)
    for i in range(1, N):
        z_mean[i] = np.sum(np.sin(i*(2*l-1)*np.pi/N)/((2*l-1)*np.pi*(np.sin((2*l-1)*np.pi/(2*N)))**2)) * f / k_H
    return z_mean

def rouse_mean2(N, T):
    f = 1./float(T)
    k_H = 3
    z_mean = np.zeros(N)
    for i in range(1, N):
        sumTerm = 0
        for m in range(1, N):
            for j in range(1, N):
                sumTerm += np.sin(i*m*np.pi/N)*np.sin(j*m*np.pi/N)/np.sin(m*np.pi/(2*N))**2
        z_mean[i] = sumTerm * f / (2*k_H*N)
    return z_mean

def rouse_var(N, T):
    k_H = 3.
    rvar = np.zeros(N)
    k = np.arange(1, N)
    for i in range(1, N):
        rvar[i] = 3.*np.sum((np.sin((i*k*np.pi)/N)/np.sin((k*np.pi)/(2*N)))**2)/(2*k_H*N)
    return rvar

def rouse_var2(N, T):
    f = 1./float(T)
    k_H = 3
    z_var = np.zeros(N)
    lam = 4*k_H*np.sin(np.arange(N)*np.pi/(2*N))**2
    omega = np.zeros([N, N])
    for i in range(N):
        for j in range(N):
            omega[i,j]=np.sqrt(2./N)*np.sin(i*j*np.pi/N)
    for i in range(1, N):
        for m in range(1, N):
            z_var[i] += 3*omega[i,m]*omega[i,m]/lam[m]
    return z_var
      
