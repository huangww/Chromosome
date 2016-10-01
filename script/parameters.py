import numpy as np

# This is a script calculating parameters for the simulation

kB = 1.38e-23 
T = 273.15 + 25
mu = 0.89
R = 5e-8
a = 2e-7
v = 4.17e-8

gamma = 6*np.pi*mu*R
print "gamma = ", gamma

F = gamma * v
print "F = ", F

f = kB*T/a
print "f = ", f

t = gamma*a*a/(kB*T)
print "t = ", t



