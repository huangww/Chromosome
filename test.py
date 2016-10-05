import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy

# use your x,y and z arrays here
x = numpy.random.rand(100)
y = numpy.random.rand(100)
z = numpy.random.rand(100)

yy, xx = numpy.meshgrid(y,x)
zz = griddata(x,y,z,xx,yy, interp='linear')
# plt.pcolor(zz)
plt.contourf(xx,yy,zz) # if you want contour plot
# plt.imshow(zz)
plt.colorbar()

plt.show()
