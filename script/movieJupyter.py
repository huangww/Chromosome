from vpython import *
import numpy as np
import os
import glob

# load data 
dataDir = 'data/'
fname = max(glob.iglob(dataDir + '*_r_*.dat'), key=os.path.getctime)
print fname

# detect parameters from file name
parameters = fname.split('_')
Nstr = [s[1:] for s in parameters if s[0]=='N']
Tstr = [s[1:] for s in parameters if s[0]=='T']
N = int(Nstr[0])
Teff = int(Tstr[0])

# set up the scene or movie
bg = vector(1., 1., 1.)
scene = canvas(width=1200, height=680, background=bg)
factor = 0.12
scale = 0.5

# faxis = frame()
# faxis.pos = [-5, -2, 0]

# load data
data = np.loadtxt(fname,comments='#')
pos = scale*data.reshape([-1,N,3])

# generate beads
bead=[]
# monomer = [279,227,123]
bead.append(sphere(pos = vec(0,0,0),
    radius = 1.5*factor,
    # radius = 2.0*factor,
    color = color.magenta))
    # frame = faxis))
bead.extend([sphere(pos = vec(pos[0,i,0],pos[0,i,1],pos[0,i,2]),
    radius = 1.0*factor,
    color = color.green)
    # frame = faxis)
    for i in range(1,N)])
   
# Connect beads with rods
# for a simple ring 
# link = [[i, (i+1)%N] for i in range(N)]
# or load topology data
link = np.loadtxt(dataDir + 'topo.dat', dtype = 'int')
# rod = [helix(pos=bead[l[0]].pos, 
rod = [cylinder(pos=bead[l[0]].pos, 
    axis = bead[l[1]].pos - bead[l[0]].pos,
    radius = 0.6 * factor, 
    color = bead[max(l[0], l[1])].color)
    # color = color.yellow,
    # thickness = 0.1*factor,
    # coils = 8)
    # frame = faxis)
    for l in link]

# calculate the center position of the scene
center = vec((np.max(pos[:,:,0]) + np.min(pos[:,:,0]))/2. + 0.6,\
	  (np.max(pos[:,:,1]) + np.min(pos[:,:,1]))/2.,\
	  (np.max(pos[:,:,2]) + np.min(pos[:,:,2]))/2.)
scene.center = center
scene.autoscale = False

# label = label(text='',height=20,color=(1.0,0,0))
#start the mainloop
print len(pos)
t=0
scene.waitfor('click')
while t<len(pos):
# while 1:
    j = t%len(pos)
    # label.text='Frame=%1.0f'%j
    for i in range(N):
	# bead[i].pos = vector(pos[j,i])
	bead[i].pos = vector(pos[j,i,0],pos[j,i,1],pos[j,i,2])
    for i,l in zip(range(len(link)), link):
	rod[i].pos=bead[l[0]].pos
	rod[i].axis=bead[l[1]].pos-bead[l[0]].pos
    # im = ImageGrab.grab()
    # im.save('tmp/tmp.%04d.png'%j)
    # scene.waitfor('click')
    t = t+1
    rate(30)

# print 'Done'
