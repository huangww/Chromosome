from visual import *
from commands import *
import numpy as np
import povexport
import os
# import glob

def GetScreenShot(FrameNumber):        
    tmp = getoutput('screencapture -R"0,42,960,457" '"$workdir"'/temp/temp.%04d.jpg' % FrameNumber)
    return

scene = display(title='Movie',width=960, height=480, background=(1.0,1.0,1.0))
factor = 0.3
scale = 0.20

beadNumber = 100
frameNumber = 2000
step = 1000

dataDir = 'lammps/'
# generate beads
bead=[]
bead.append(sphere(pos=(0,0,0),radius=2.0*factor,color=(1.0,0.0,0.0)))
for i in range(1,beadNumber):
    bead.append(sphere(pos=(i*scale,0,0),radius=1.0*factor,color=(0.2,0.8,0.2)))

# Connect beads 
link = np.loadtxt(dataDir + 'topol.dat', dtype = 'int')
rod = []
for l in link:
    if (l[1]>l[0]):
        rod.append(cylinder(pos=bead[l[0]].pos,axis=bead[l[1]].pos-bead[l[0]].pos,radius=0.5*factor,color=bead[l[1]].color))
    else:
        rod.append(cylinder(pos=bead[l[0]].pos,axis=bead[l[1]].pos-bead[l[0]].pos,radius=0.5*factor,color=bead[l[0]].color))


#start the mainloop
scene.waitfor('click keydown')
t=0
while t<=frameNumber:
    datafile = dataDir + 'dump' + str(t*step) + '.dat'
    beadPosition = np.loadtxt(datafile, skiprows = 2, usecols = (1, 2, 3))
    for i in range(beadNumber):
        bead[i].pos= vector(beadPosition[i])
    i=0
    for l in link:
        rod[i].pos=bead[l[0]].pos
        rod[i].axis=bead[l[1]].pos-bead[l[0]].pos
        i = i+1
    rate(50)
    #GetScreenShot(t)
    t = t+1
    # scene.waitfor('click keydown')
    # s = scene.kb.getkey()
    # if s == 's':
        # povexport.export(display = scene)

print 'Done'
	
