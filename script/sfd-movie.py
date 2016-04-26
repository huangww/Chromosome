import numpy as np
import matplotlib.pyplot as plt
from commands import *
from visual import *
import povexport
import os
import glob

def GetScreenShot(FrameNumber):        
    tmp = getoutput('screencapture -R"0,42,800,300" fig/temp.%04d.png' % FrameNumber)
    return


# load data 
dataDir = 'data/'
fname = max(glob.iglob(dataDir + 'par*.dat'), key=os.path.getctime)
# fname = dataDir + 'par_N'+str(N)+'_T' + str(Teff) + '.dat'
print fname

# detect parameters from file name
parameters = fname.split('_')
Nstr = [s[1:] for s in parameters if s[0]=='N']
Tstr = [s[1:-4] for s in parameters if s[0]=='T']
N = int(Nstr[0])
Teff = int(Tstr[0])
print 'N = ', N
print 'Teff = ', Teff

parPos = np.loadtxt(fname, comments = '#')
# state = data.reshape([-1, N, 2])

# set up the scene or movie
scene = display(width=800, height=600, background=(1.0,1.0,1.0))
scene.fullscreen = True
scene.autoscale = False
scale = 0.2

faxis = frame()
faxis.pos = (-10, 2, 0)
axis = curve(x = arange(N+1), y = 0, z=0, color=color.blue, frame = faxis)
ticks = [curve(pos=[(i, 0, 0), (i, 0.2, 0)], 
    color = color.blue,
    frame = faxis) 
    for i in range(N+1)]
particles = [sphere(pos=vector(parPos[0,i+1], 0, 0), 
    radius = 1.0*scale,  
    color=(0.2, 0.8, 0.2), 
    frame=faxis)
    for i in range(N/2)]

t=0
while t < len(parPos):
    j = t%len(parPos)
    for par in particles:
        par.pos = vector(parPos[j,particles.index(par)+1],0,0)
    if t == 0:
        s = scene.kb.getkey()
    GetScreenShot(t)
    print t
    t = t+1
    rate(30)

