import numpy as np
import matplotlib.pyplot as plt
from theory import *

fig = plt.figure(0)
ax = fig.add_subplot(111)

N = 100
T = 1000


index = np.arange(0,101)
var = np.zeros(len(index))
varbench = np.zeros(len(index))
for i in index:
    varbench[i] = i*(N-i)/float(N)

ax.plot(index, varbench, 'ro-')
for T in [10,20,100,10000]:
    for i in index:
        var[i] = r_var(i, N, T)
    ax.plot(index, var, label = 'T='+str(T))

ax.legend()
plt.show()

