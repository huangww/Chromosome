from vpython import *
import numpy as np

scene = canvas(width=800, height=480, background=vector(1,1,1))
L = 10
lattice = cylinder(pos=vector(0,0,0), axis=vector(L, 0, 0), radius=0.1, color=color.cyan, opacity=0.2)
ticks = [cylinder(pos=vector(i, -0.1, 0), axis=vector(0, 0.2, 0), radius = 0.02, color = color.blue) for i in range(L+1)]
wallLeft = box(pos=vector(-0.05,0,0), length=0.05, height=1, width=1, color=color.yellow, opacity=0.8)
wallRight = box(pos=vector(L+0.05,0,0), length=0.05, height=1, width=1, color=color.yellow, opacity=0.8)

# loop = ring(pos=vector(5,0,0), axis=vector(0,1,0), radius=5, thickness=0.1)
pointerRight = arrow(pos=vector(4,1,0),axis=vector(3,0,0), shaftwidth=0.2, color=color.magenta, opacity=0.3)
pointerLeft = arrow(pos=vector(7,-1,0),axis=vector(-3,0,0), shaftwidth=0.2, color=color.magenta, opacity=0.3)

N = L / 2
initPos = [0, 2, 5, 6, L-1]
particle = [sphere(pos=vector(initPos[i], 0, 0), radius = 0.2,  color=color.green, opacity=0.8) for i in range(N)]
particle[0].color = color.red


# calculate the postion of each particle 
t = 0
tEnd = 100
neighbor = 1
posFrame = np.zeros([tEnd+1, N])
posFrame[0] = initPos
while t < tEnd:
    t = t + 1
    posFrame[t] = posFrame[t-1]
    i = t % (2*L)
    i = 2*L-i if i > L else i
    i0 = (t-1) % (2*L)
    i0 = 2*L-i0 if i0 > L else i0
    step = i - i0
    posFrame[t, 0] = i
    if posFrame[t, neighbor] == posFrame[t, 0]:
        posFrame[t, neighbor] -= step 
        neighbor += step
        neighbor = N-1 if neighbor==N else neighbor
        neighbor = 1 if neighbor==0 else neighbor

t=1
scene.center = vector(5., 0, 0)
scene.waitfor('click')
while t < tEnd+1:
    for i in range(N):
        particle[i].pos = vector(posFrame[t, i], 0, 0)
    t = t + 1
    rate(2)

