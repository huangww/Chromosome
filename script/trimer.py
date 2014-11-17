import numpy as np
import numpy.linalg as la
import matplotlib.pylab as plt

dataDir = 'data/trimer/'
beadNumber = 3
data = np.loadtxt(dataDir + 'trimer.dat')
pos = data.reshape([3,beadNumber,-1],order='F')

vec1 = pos[:,0]-pos[:,1]
vec2 = pos[:,2]-pos[:,1]
cosang = np.zeros(len(vec1[0]))
#sinang = np.zeros(len(vec1[0]))
for i in range(len(vec1[0])):
    cosang[i] = np.dot(vec1[:,i],vec2[:,i])
    #sinang[i] = la.norm(np.cross(vec1[:,i],vec2[:,i]))
#theta = np.arctan2(sinang, cosang)
plt.subplot(121)
plt.hist(cosang,bins=50,normed=True)
x = np.linspace(0,np.pi,100)
plt.plot(np.cos(x),0.5*np.ones(len(x)),linewidth=3,label='Bead Spring')
plt.plot(np.cos(x),0.523*np.sqrt(1-0.25*np.cos(x)*np.cos(x)),linewidth=3,label='Bead Rod')
plt.ylim([0.4,0.6])
plt.xlabel(r'$cos\theta$')
plt.ylabel(r'$p(cos\theta)$')
plt.title('Without Pseudo Force')
plt.legend()

data = np.loadtxt(dataDir + 'trimer_pseudo.dat')
pos = data.reshape([3,beadNumber,-1],order='F')
vec1 = pos[:,0]-pos[:,1]
vec2 = pos[:,2]-pos[:,1]
cosang = np.zeros(len(vec1[0]))
for i in range(len(vec1[0])):
    cosang[i] = np.dot(vec1[:,i],vec2[:,i])

plt.subplot(122) 
plt.hist(cosang,bins=50,normed=True)
x = np.linspace(0,np.pi,100)
plt.plot(np.cos(x),0.5*np.ones(len(x)),linewidth=3,label='Bead Spring')
plt.plot(np.cos(x),0.523*np.sqrt(1-0.25*np.cos(x)*np.cos(x)),linewidth=3,label='Bead Rod')
plt.ylim([0.4,0.6])
plt.xlabel(r'$cos\theta$')
plt.ylabel(r'$p(cos\theta)$')
plt.title('With Pseudo Force')
plt.legend()

plt.show()

