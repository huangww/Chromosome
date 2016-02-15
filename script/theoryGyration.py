import numpy as np

def e_z_mean(i, N, T):
    mu = (N+1)/2.0
    x = (mu - i)/float(T)
    mean_e_z = 1.0/np.tanh(x) - 1.0/x
    return mean_e_z

def e_z_var(i, N, T):
    mu = (N+1)/2.0
    x = (mu - i)/float(T)
    var_e_z = 1/x**2 - 1/np.sinh(x)**2 
    return var_e_z

def Gyration(N, T):
    index = np.arange(N)
    e_z = e_z_mean(index, N, T)
    var_e_z = e_z_var(index, N, T)
    second_moment_e_z = var_e_z + e_z*e_z
    var_e_xy = 1 - second_moment_e_z
    sum_var_e_z = var_e_z.sum()
    sum_var_e_xy = var_e_xy.sum()
    z_mean = np.zeros(len(index))
    z_covar = np.zeros([len(index),len(index)])
    xy_covar = np.zeros([len(index),len(index)])
    for i in index:
        z_mean[i] = (e_z[0:i]).sum()
        for j in range(i,len(index)):
            z_covar[i,j] = (var_e_z[:i]).sum()* \
                    (var_e_z[j:]).sum()/sum_var_e_z
            xy_covar[i,j] = sum(var_e_xy[:i])* \
                    sum(var_e_xy[j:])/sum_var_e_xy
            z_covar[j,i] = z_covar[i,j]
            xy_covar[j,i] = xy_covar[i,j]
    z_mean_ij = np.outer(z_mean,z_mean)
    lam_z = (z_covar.trace()+(z_mean*z_mean).sum())/N \
            - (z_covar+z_mean_ij).sum()/(N*N)
    lam_xy = xy_covar.trace()/N - xy_covar.sum()/(N*N)
    return lam_z + lam_xy

def GyrationTensorNumerical(beadPos):
    N = len(beadPos[0])
    cm = beadPos.mean(axis = 1)
    centredPos = beadPos
    for i in range(N):
        centredPos[:,i,:] = beadPos[:,i,:] - cm

    Sxx = centredPos[:,:,0]*centredPos[:,:,0]
    Sxx = Sxx.sum(axis=1).mean()/N
    Sxy = centredPos[:,:,0]*centredPos[:,:,1]
    Sxy = Sxy.sum(axis=1).mean()/N
    Sxz = centredPos[:,:,0]*centredPos[:,:,2]
    Sxz = Sxz.sum(axis=1).mean()/N
    Syx = centredPos[:,:,1]*centredPos[:,:,0]
    Syx = Syx.sum(axis=1).mean()/N
    Syy = centredPos[:,:,1]*centredPos[:,:,1]
    Syy = Syy.sum(axis=1).mean()/N
    Syz = centredPos[:,:,1]*centredPos[:,:,2]
    Syz = Syz.sum(axis=1).mean()/N
    Szx = centredPos[:,:,2]*centredPos[:,:,0]
    Szx = Szx.sum(axis=1).mean()/N
    Szy = centredPos[:,:,2]*centredPos[:,:,1]
    Szy = Szy.sum(axis=1).mean()/N
    Szz = centredPos[:,:,2]*centredPos[:,:,2]
    Szz = Szz.sum(axis=1).mean()/N

    S = np.array([[Sxx,Sxy,Sxz],[Syx,Syy,Syz],[Szx,Szy,Szz]])
    return S
    
def GyrationTensor(N, T):
    index = np.arange(N)
    e_z = e_z_mean(index, N, T)
    var_e_z = e_z_var(index, N, T)
    second_moment_e_z = var_e_z + e_z*e_z
    var_e_xy = 1 - second_moment_e_z
    sum_var_e_z = var_e_z.sum()
    sum_var_e_xy = var_e_xy.sum()
    z_mean = np.zeros(len(index))
    z_covar = np.zeros([len(index),len(index)])
    xy_covar = np.zeros([len(index),len(index)])
    for i in index:
        z_mean[i] = (e_z[0:i]).sum()
        for j in range(i,len(index)):
            z_covar[i,j] = (var_e_z[:i]).sum()* \
                    (var_e_z[j:]).sum()/sum_var_e_z
            xy_covar[i,j] = sum(var_e_xy[:i])* \
                    sum(var_e_xy[j:])/sum_var_e_xy
            z_covar[j,i] = z_covar[i,j]
            xy_covar[j,i] = xy_covar[i,j]
    z_mean_ij = np.outer(z_mean,z_mean)
    lam_z = (z_covar.trace()+(z_mean*z_mean).sum())/N \
            - (z_covar+z_mean_ij).sum()/(N*N)
    lam_xy = xy_covar.trace()/N - xy_covar.sum()/(N*N)
    S = np.zeros([3,3])
    S[0,0] = lam_z
    S[1,1] = lam_xy/2.
    S[2,2] = lam_xy/2.
    return S
    
def Asphericity(S):
    Q = S - S.trace()/3*np.identity(3)
    Qs = np.dot(Q,Q)
    # kappa = (S[0,0]*S[0,0]+S[1,1]*S[1,1]+S[2,2]*S[2,2]) \
    #         /S.trace()**2 * 1.5 - 0.5
    kappa = 3*Qs.trace()/(2*S.trace()**2)
    sigma = (4*np.linalg.det(Q))/((2.*Qs.trace()/3.)**(3./2))
    return kappa, sigma

def AsphericityArray(beadPos):
    kappa = np.zeros(len(beadPos))
    sigma = np.zeros(len(beadPos))
    count = 0
    for frame in beadPos:
        cm = frame.mean(axis = 0)
        centredPos = frame - cm
        Sxx = centredPos[:,0]*centredPos[:,0]
        Sxx = Sxx.mean()
        Sxy = centredPos[:,0]*centredPos[:,1]
        Sxy = Sxy.mean()
        Sxz = centredPos[:,0]*centredPos[:,2]
        Sxz = Sxz.mean()
        Syx = centredPos[:,1]*centredPos[:,0]
        Syx = Syx.mean()
        Syy = centredPos[:,1]*centredPos[:,1]
        Syy = Syy.mean()
        Syz = centredPos[:,1]*centredPos[:,2]
        Syz = Syz.mean()
        Szx = centredPos[:,2]*centredPos[:,0]
        Szx = Szx.mean()
        Szy = centredPos[:,2]*centredPos[:,1]
        Szy = Szy.mean()
        Szz = centredPos[:,2]*centredPos[:,2]
        Szz = Szz.mean()
        S = np.array([[Sxx,Sxy,Sxz],[Syx,Syy,Syz],[Szx,Szy,Szz]])
        kappa[count], sigma[count] = Asphericity(S)
        count = count + 1
    return kappa, sigma

