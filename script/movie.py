from visual import *
from commands import *
import numpy as np
import povexport
import os
import glob

def GetScreenShot(FrameNumber):        
    tmp = getoutput('screencapture -R"0,42,960,457" '"$workdir"'/temp/temp.%04d.jpg' % FrameNumber)
    return

#load the newest data
dataDir = 'data/'
fileName = max(glob.iglob(dataDir + '*.dat'), key=os.path.getctime)
data = np.loadtxt(fileName)
#get parameters from the fileName
parameters = fileName.split('_')
for s in parameters:
    if s[0] == 'N':
        beadNumber = int(s[1:])

scene = display(title='Title',width=960, height=480, background=(1.0,1.0,1.0))
factor = 0.3
scale = 0.20
# generate beads
bead=[]
bead.append(sphere(pos=(0,0,0),radius=2.0*factor,color=(1.0,0.0,0.0)))
for i in range(1,beadNumber):
    bead.append(sphere(pos=(i*scale,0,0),radius=1.0*factor,color=(0.8,0.2,0.8)))
    #i0 = i%monomer.sum()
    #if (i0>0 and i0<=monomer[0]):
        #bead.append(sphere(pos=(i*scale,0,0),radius=1.0*factor,color=(0.8,0.2,0.8)))
    #elif (i0>0 and i0<=monomer[0]+monomer[1]):
        #bead.append(sphere(pos=(i*scale,0,0),radius=1.0*factor,color=(0.2,0.8,0.2)))
    #else: 
    #    bead.append(sphere(pos=(i*scale,0,0),radius=1.0*factor,color=(0.2,0.8,0.8)))

# Connect beads with rods
# for a simple ring 
#link = np.zeros((beadNumber,2),dtype=int)
#for i in range(beadNumber):
    #link[i,0]=i
#    link[i,1]=(i+1)%beadNumber

link = np.loadtxt(dataDir + 'topol.dat', dtype = 'int')
rod = []
for l in link:
    if (l[1]>l[0]):
        rod.append(cylinder(pos=bead[l[0]].pos,axis=bead[l[1]].pos-bead[l[0]].pos,radius=0.5*factor,color=bead[l[1]].color))
    else:
        rod.append(cylinder(pos=bead[l[0]].pos,axis=bead[l[1]].pos-bead[l[0]].pos,radius=0.5*factor,color=bead[l[0]].color))

beadPosition = data.reshape([3,beadNumber,-1],order='F')
#temp = np.reshape(data,(2,beadNumber,-1),order='F')
#beadPosition = np.zeros((3,beadNumber,len(temp[0,0])))
#beadPosition[0,:,:]=temp[0,:,:]
#beadPosition[1,:,:]=temp[1,:,:]

#calculate the center position of the system
center = [np.max(beadPosition[0,:,:])/2.+np.min(beadPosition[0,:,:])/2.,\
	  np.max(beadPosition[1,:,:])/2.+np.min(beadPosition[1,:,:])/2.,\
	  np.max(beadPosition[2,:,:])/2.+np.min(beadPosition[2,:,:])/2.]
scene.center=center
#scene.autoscale=False
#label = label(pos=[center[0],np.min(beadPosition[1,:,:]),center[2]],text='',height=20,color=(1.0,0,0))

#start the mainloop
scene.waitfor('click keydown')
t=0
#while t<len(beadPosition[0,0]):
while 1:
    #label.text='Frame=%1.0f'%t
    j = t%len(beadPosition[0,0])
    for i in range(beadNumber):
	bead[i].pos= vector(beadPosition[:,i,j])
    i=0
    for l in link:
	rod[i].pos=bead[l[0]].pos
	rod[i].axis=bead[l[1]].pos-bead[l[0]].pos
	i = i+1
    rate(50)
    #GetScreenShot(t)
    t = t+1
    scene.waitfor('click keydown')
    s = scene.kb.getkey()
    if s == 's':
        povexport.export(display = scene)

print 'Done'
	
