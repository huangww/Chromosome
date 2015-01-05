import numpy as np
import matplotlib.pyplot as plt

theta = np.linspace(0, np.pi/3., 100)

rho = -1/np.cos(theta + 2*np.pi/3)


plt.plot(theta, rho)
plt.xlim([0,np.pi/3])
plt.ylim([0,2])


plt.show()
