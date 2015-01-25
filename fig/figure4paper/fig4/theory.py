import numpy as np

def integral_xy_var(i, N, T):
    mu = (N+1)/2.0
    mean_z = 1.0/np.tanh((mu-i)/T)-T/(mu-i)
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
    z_var = T*((1.0/np.tanh((n-mu)/T)-1.0/((n-mu)/T))-(1.0/np.tanh((m-mu)/T)-1.0/((m-mu)/T)))
    return z_var

def var_z(N, T):
    var_z = np.zeros(N)
    for i in range(N):
        v0s = z_var(0,i,N,T)
        vst = z_var(i,N,N,T)
        v0t = v0s + vst
        var_z[i] = 2.0*v0s*vst/v0t
    return var_z

def var_xy(N, T):
    var_xy = np.zeros(N)
    for i in range(N):
        v0s = xy_var(0,i,N,T)
        vst = xy_var(i,N,N,T)
        v0t = v0s + vst
        var_xy[i] = 2.0*v0s*vst/v0t
    return var_xy

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

