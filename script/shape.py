from visual import *

scene = display(width = 960, height = 480, background = (1,1,1))
ob1 = sphere(pos=(0,0,0), radius = 1.2, color = color.green)

ob2 = ellipsoid(pos=(6,0,0), axis=(1,0,0), \
                 length=4, height=2, width=2, color = color.green)

ob3 = rod = cylinder(pos=(10,0,0), axis=(4,0,0), radius=0.2, color = color.green)
