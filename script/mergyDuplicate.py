import numpy as np


def Mergy(data):
    keyArray = np.unique(data[:,0])
    valArray = np.zeros(len(keyArray))
    for key,valIndex in zip(keyArray, range(len(valArray))):
        val = 0
        count = 0
        for i in range(len(data)):
            if data[i,0] == key:
                val = val + data[i,1]
                count = count +1
        val = val / count
        valArray[valIndex] = val
    return np.transpose([keyArray, valArray])

def main():
    data = np.loadtxt('data/tauTemp3D.dat')
    data = Mergy(data)
    np.savetxt('data/tauForce3D.dat', data)

if __name__ == "__main__":
    main()
