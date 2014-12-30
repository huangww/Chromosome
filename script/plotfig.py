import matplotlib.pyplot as plt
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

# initialize the figure
fig = plt.figure(0)
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)
font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 12 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text',usetex = True)

beadNumber = 500
# frameNumber = 5000
dimension = 3
index = np.linspace(0,beadNumber,beadNumber)
tempFrames = 0
step = 1

# for T in ['1', '5','10','5000']:
for T in ['1']:
    # load data
    dataDir = 'data/'
    # fileName = dataDir + 'r_N' + str(beadNumber) + \
            # '_T'+str(T) + '_5489.dat'
    fileName = dataDir + 'MD_N' + str(beadNumber) + \
            '_T'+str(T) + '_sample2.dat'
    data = np.loadtxt(fileName,comments='#')
    beadPos = data.reshape([-1,beadNumber,3])
    
    # calculate mean and varriance from data
    posMean = beadPos.mean(axis = 0)
    posVars = beadPos.var(axis = 0)

    # plot figure
    N = beadNumber
    line, = ax1.plot(index, z_mean(index,N,int(T)))
    ax1.plot(index, posMean[:,0],line.get_color()+'o')
    line, = ax2.plot(index, z_var(index,N,int(T)))
    ax2.plot(posVars[:,0],line.get_color()+'o')
    ax3.plot(index, posMean[:,1],line.get_color()+'o')
    ax3.plot(index, posMean[:,2],line.get_color()+'*')
    var_xy = np.zeros(len(index))
    for i in range(N):
        var_xy[i] = xy_var(i,N,int(T))
    line, = ax4.plot(index, var_xy/2.0)
    ax4.plot(index, posVars[:,1],line.get_color()+'o')
    ax4.plot(index, posVars[:,2],line.get_color()+'*')

ax3.set_ylim([-5,5])

plt.show()

