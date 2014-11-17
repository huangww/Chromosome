import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob as glob

def ExportData():
    filelist = glob.glob('*.trcT')
    for pair in range(len(filelist)/2):
        x1 = []
        y1 = []
        inputfile = open(filelist[2*pair], 'r')
        for line in inputfile:
            line = line.strip()
            if line[0] != '%':
                columns = line.split()
                x1.append(float(columns[2]))
                y1.append(float(columns[3]))
        x2 = []
        y2 = []
        inputfile = open(filelist[2*pair+1], 'r')
        for line in inputfile:
            line = line.strip()
            if line[0] != '%':
                columns = line.split()
                x2.append(float(columns[2]))
                y2.append(float(columns[3]))
        np.savetxt(str(pair+1)+'.txt',(x1,y1,x2,y2))

        fig = plt.figure()
        plt.subplot(121)
        x1 = np.array(x1)
        y1 = np.array(y1)
        x2 = np.array(x2)
        y2 = np.array(y2)
        plt.plot(np.sqrt(x1*x1+y1*y1))
        plt.plot(np.sqrt(x2*x2+y2*y2))
        plt.plot(np.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)))
        plt.legend(['$r_1$','$r_2$','$r_{12}$'])
        plt.subplot(122)
        for i in range(len(x1)-1):
            p1 = plt.plot([x1[i],x1[i+1]], [y1[i],y1[i+1]], '-o', linewidth=2, markersize=8, color=plt.cm.rainbow(float(i)/float(np.size(x1))),alpha=0.5)
            p2 = plt.plot([x2[i],x2[i+1]], [y2[i],y2[i+1]], '-v', linewidth=2, markersize=8,color=plt.cm.rainbow(float(i)/float(np.size(x1))),alpha=1.0)
        lg = plt.legend(["ADE_1","ADE_2"])
        ax2 = fig.add_axes([0.9,0.1,0.05,0.8])
        cb = mpl.colorbar.ColorbarBase(ax2,cmap=plt.cm.rainbow )
        fig.savefig(str(pair+1)+'.png')
        plt.close(fig)
    return





