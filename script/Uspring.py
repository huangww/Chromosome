import numpy as np
import matplotlib.pyplot as plt

r0 = 1.1
sigma = 0.9
k = 1.0
eps = 7.2

r = np.linspace(0.8, 1.099999, 100)
u = -0.5 * k * r0 * r0 * np.log(1-(r/r0)**2) + \
        4 * eps * ((sigma/r)**12 - (sigma/r)**6) + eps

fig = plt.plot(r, u)

plt.show()
