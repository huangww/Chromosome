import matplotlib.pyplot as plt

font = {'family' : 'sans-serif',
        'serif'  : 'Helvetica',
        'weight' : 'normal',
        'size'   : 20 }
plt.rc('lines', lw=2)
plt.rc('font', **font)
plt.rc('text', usetex=True)


Dir = './'
fileName = Dir + 'ring.png'
image = plt.imread(fileName)

plt.imshow(image)
plt.text(100,80, r'$x$')
plt.arrow(100, 300, 0, -200, head_width=5, \
        head_length=8, lw = 2, fc='k', ec='k')
plt.text(30,430, r'$y$')
plt.arrow(100, 300, -50, 100, head_width=5, \
        head_length=8, lw = 2, fc='k', ec='k')
plt.text(600,320, r'$z$')
plt.arrow(100, 300, 500, 0, head_width=5, \
        head_length=8, lw = 2, fc='k', ec='k')

plt.text(315,155, r'$\mathbf{r}_i$')
plt.text(330,220, r'$\mathbf{e}_i$')
plt.text(365,175, r'$\mathbf{r}_{i+1}$')
plt.arrow(319.5, 175, 53, 17, head_width=5, \
        head_length=8, lw = 2, fc='k', ec='k')

plt.text(370,110, r'$\mathbf{F}$')
plt.arrow(250, 100, 100, 0, head_width=5, \
        head_length=8, lw = 2, fc='k', ec='k')

plt.text(395,285, r'$a$')
plt.plot([362.081, 430],[270.468,255],'k-o')

plt.axis('off')
plt.savefig('figure1c.pdf')
plt.show()
