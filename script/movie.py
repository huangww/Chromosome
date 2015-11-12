from visual import *
from commands import *
import numpy as np
import povexport
import os
import glob

def GetScreenShot(FrameNumber):        
    tmp = getoutput('screencapture -R"0,42,960,457" data/temp.%04d.png' % FrameNumber)
    return

# load data 
dataDir = 'data/'
# fname = max(glob.iglob(dataDir + 'r*.dat'), key=os.path.getctime)
fname = dataDir + 'r_N100_T1_1.dat'
print fname

# detect parameters from file name
parameters = fname.split('_')
Nstr = [s[1:] for s in parameters if s[0]=='N']
Tstr = [s[1:] for s in parameters if s[0]=='T']
N = int(Nstr[0])
Teff = int(Tstr[0])

# set up the scene or movie
scene = display(width=1200, height=600, background=(1.0,1.0,1.0))
scene.fullscreen = True
scene.autoscale = False
# calculate the center position of the scene
# center = [(np.max(pos[:,:,0]) + np.min(pos[:,:,0]))/2.,\
# 	  (np.max(pos[:,:,1]) + np.min(pos[:,:,1]))/2.,\
# 	  (np.max(pos[:,:,2]) + np.min(pos[:,:,2]))/2.]
# scene.center = center
# print center

factor = 0.15
scale = 0.4

faxis = frame()
faxis.pos = [-11, 0, 0]

# load data
data = np.loadtxt(fname,comments='#')
pos = scale*data.reshape([-1,N,3])


# generate beads
bead=[]
# monomer = [279,227,123]
bead.append(sphere(pos = (0,0,0),
    radius = 2.0*factor,
    color = color.red,
    frame = faxis))
bead.extend([sphere(pos = (i,0,0),
    radius = 1.0*factor,
    color = color.green,
    frame = faxis)
    for i in range(1,N)])
   
# Connect beads with rods
# for a simple ring 
link = [[i, (i+1)%N] for i in range(N)]
# or load topology data
# link = np.loadtxt(dataDir + 'topo.dat', dtype = 'int')
rod = [cylinder(pos=bead[l[0]].pos, 
    axis = bead[l[1]].pos - bead[l[0]].pos,
    radius = 0.5 * factor, 
    color = bead[max(l[0], l[1])].color,
    frame = faxis)
    for l in link]


# label = label(text='',height=20,color=(1.0,0,0))
#start the mainloop
t=0
# while t<len(beadPosition):
while 1:
    j = t%len(pos)
    # label.text='Frame=%1.0f'%j
    for i in range(N):
	bead[i].pos = vector(pos[j,i])
    for i,l in zip(range(N), link):
	rod[i].pos=bead[l[0]].pos
	rod[i].axis=bead[l[1]].pos-bead[l[0]].pos
    rate(30)
    if t == 0:
        scene.waitfor('click keydown')
    # GetScreenShot(t)
    # if t>100:
    #     s = scene.kb.getkey()
    # if s == 's':
        # povexport.export(display = scene)
    t = t+1

print 'Done'
