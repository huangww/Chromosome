import numpy as np
import numpy.linalg as la
import matplotlib.pylab as plt

dataDir = 'data/trimer/'
beadNumber = 3

def getCosTheta(data):
    pos = data.reshape([-1,beadNumber,3])
    vec1 = pos[:,0]-pos[:,1]
    vec2 = pos[:,2]-pos[:,1]
    cosang = np.zeros(len(vec1))
    #sinang = np.zeros(len(vec1[0]))
    for t in range(len(vec1)):
        cosang[t] = np.dot(vec1[t],vec2[t])
        if cosang[t]>1.:
            cosang[t] = 1.
        if cosang[t]<-1.:
            cosang[t] = -1.
        #sinang[i] = la.norm(np.cross(vec1[:,i],vec2[:,i]))
    #theta = np.arctan2(sinang, cosang)
    return cosang

def plotfig():
    fig = plt.figure(0)
    data = np.loadtxt(dataDir + 'trimer.dat')
    cosang = getCosTheta(data)
    ax = fig.add_subplot(121)
    ax.hist(cosang,bins=30,normed=True)
    x = np.linspace(0,np.pi,100)
    ax.plot(np.cos(x),0.5*np.ones(len(x)),linewidth=3,label='Uniform')
    ax.plot(np.cos(x),0.523*np.sqrt(1-0.25*np.cos(x)*np.cos(x)),linewidth=3,label='Non-uniform')
    ax.set_xlim([-1,1])
    ax.set_ylim([0.4,0.6])
    ax.set_xlabel(r'$\cos\theta$')
    ax.set_ylabel(r'$p(\cos\theta)$')
    ax.set_title('Without Pseudo Force')
    # axinset = fig.add_axes([0.2,0.6,0.25,0.25])
    # image = plt.imread('fig/sketch/trimer.png')
    # axinset.imshow(image)
    # axinset.axis('off')
    # ax.legend()

    ax = fig.add_subplot(122)
    data = np.loadtxt(dataDir + 'trimer_pseudo.dat')
    cosang = getCosTheta(data)
    ax.hist(cosang,bins=30,normed=True)
    x = np.linspace(0,np.pi,100)
    ax.plot(np.cos(x),0.5*np.ones(len(x)),linewidth=3,label='Uniform')
    ax.plot(np.cos(x),0.523*np.sqrt(1-0.25*np.cos(x)*np.cos(x)),linewidth=3,label='Non-uniform')
    ax.set_xlim([-1,1])
    ax.set_ylim([0.4,0.6])
    ax.set_xlabel(r'$\cos\theta$')
    ax.set_ylabel(r'$p(\cos\theta)$')
    ax.set_title('With Pseudo Force')
    ax.legend()
    plt.show()


def main():
    plotfig()

if __name__ == "__main__":
    main()
