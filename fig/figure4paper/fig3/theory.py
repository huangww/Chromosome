import numpy as np

def mean_e_z(i, N, T):
    mu = (N+1)/2.0
    x = (mu - i)/float(T)
    mean_e_z = 1.0/np.tanh(x) - 1.0/x
    return mean_e_z

def var_e_z(i, N, T):
    mu = (N+1)/2.0
    x = (mu - i)/float(T)
    var_e_z = 1/x**2 - 1/np.sinh(x)**2 
    return var_e_z
    
def z_mean(index, N, T):
    e_z_i = mean_e_z(index, N, T)
    z_mean = np.zeros(len(index))
    for i in index:
        z_mean[i] = sum(e_z_i[0:i])
    return z_mean

def z_var(index, N, T):
    var_e_z_i = var_e_z(index, N, T)
    z_var = np.zeros(len(index))
    for i in index:
        z_var[i] = sum(var_e_z_i[:i])* \
                sum(var_e_z_i[i:])/sum(var_e_z_i)
    return z_var

def var_e_xy(index, N, T):
    meanez = mean_e_z(index, N, T)
    varez = var_e_z(index, N, T)
    second_moment_e_z = varez + meanez*meanez
    var_e_xy = 1 - second_moment_e_z
    return var_e_xy

def xy_var(index, N, T):
    var_e_xy_i = var_e_xy(index, N, T)
    xy_var = np.zeros(len(index))
    for i in index:
        xy_var[i] = sum(var_e_xy_i[:i])* \
                sum(var_e_xy_i[i:])/sum(var_e_xy_i)
    return xy_var

def r_var(index, N, T):
    r_var = xy_var(index,N,T) + z_var(index,N,T)
    return r_var

# def xy_covar(i, j, N, T):
#     s = min(i,j)
#     l = max(i,j)
#     xy_covar = xy_raw_var(0,s,N,T)*xy_raw_var(l,N,N,T)/xy_raw_var(0,N,N,T)
#     return xy_covar
#
# def z_covar(i, j, N, T):
#     s = min(i,j)
#     l = max(i,j)
#     z_covar = (z_raw_var(s,N,T)-z_raw_var(0,N,T))*(z_raw_var(N,N,T)-z_raw_var(l,N,T))/(z_raw_var(N,N,T)-z_raw_var(0,N,T))
#     return z_covar
#  
# def r_covar(i, j, N, T):
#     r_covar = xy_covar(i,j,N,T) + z_covar(i,j,N,T)
#     return r_covar
#
# def Gyration(N, T):
#     rgs = 0
#     for i in range(N):
#         rgs += r_var(i,N,T) + z_mean(i,N,T)**2
#     rgs = rgs/N
#     for i in range(N):
#         for j in range(N):
#             rgs -= (r_covar(i,j,N,T) + z_mean(i,N,T)*z_mean(j,N,T))/(N*N)
#     return rgs
#
