from __future__ import division
from sympy import *
import numpy as np
x, x1, x2, z, z1, z2, p, q, p1, p2, B1, B2, B, a, b, E = symbols('x x1 x2 z z1 z2 p q p1 p2 B1 B2 B a b E')
f, g, h = symbols('f g h', cls=Function)

L = 4
f = ((a-b*z1)*z2**x2 + (a-b*z2)*z1**x2) * (z1**(-x1)*(z2**(-L-1)-z2**(-L)) + z2**(-x1)*(z1**(-L-1)-z**(-L)))
g = (a/b)**(x2-x1-L)*(z2**(-x2)*(1/z1-1)+z1**(-x2)*(1/z2-1))*((a-b*z2)*z2**L*z1**x1 +(a-b*z1)*z1**L*z2**x1)
h = (g - f).subs([(a, 1), (b, 2)])
h1 = h.subs([(x1,1),(x2,2)])
h2 = h.subs([(x1,L-1),(x2,L)])
zr = solve([h1, h2], [z1, z2])
for val in zr:
    ev = -2*a - 2*b + a/val[0] + b*val[0] + a/val[1] + b*val[1] 
    print ev.subs([(a,1), (b, 2)]).evalf()

# f = z**x + B*(a/b)**x*z**(-x)
# g = a*(f.subs(x, 2) + a/b*f.subs(x, 0)) - b*(f.subs(x, 2) + a/b*f.subs(x, 1))
# h = a*(f.subs(x, L) + a/b*f.subs(x, L-1)) - b*(f.subs(x, L+1) + a/b*f.subs(x, L-1))

# f = z**x + B*(a/b)**x*z**(-x)
# f1 = f.subs([(z, z1), (B, B1)])
# f2 = f.subs([(z, z2), (B, B2)])
# E = a/z1 + a/z2 + b*z1 + b*z2 -2*(a-b)
# g = E*(f1.subs(x, 1)*f2.subs(x, 2) + f2.subs(x, 1)*f1.subs(x,2)) + \
#         a*(f1.subs(x, 1)*f2.subs(x, 2) + f2.subs(x, 1)*f1.subs(x,2)) \
#     - b*(f1.subs(x, 1)*f2.subs(x, 3) + f2.subs(x, 1)*f1.subs(x, 3))
# h = E*(f1.subs(x, L-1)*f2.subs(x, L) + f2.subs(x, L-1)*f1.subs(x,L)) + \
#         b*(f1.subs(x, L-1)*f2.subs(x, L) + f2.subs(x, L-1)*f1.subs(x,L)) \
#     - a*(f1.subs(x, L-2)*f2.subs(x, L) + f2.subs(x, L-2)*f1.subs(x, L))
# g = g.subs([(a, 1), (b, 2)])
# h = h.subs([(a, 1), (b, 2)])

# B0 = solve(g, B)[0]
# h = h.subs(B, B0)
# theta =  solve(h.subs([(a,1), (b,2)]), p)
# print theta
# for val in theta:
#     ev = -a-b+2*sqrt(a*b)*cos(val)
#     print ev.subs([(a,1), (b, 2)]).evalf()

# f = q*(x**L + (x-q)/(x-1)*(q/x)**(L-1)) - (x**(L+1) + (x-q)/(x-1)*(q/x)**L)
# roots = solve(f.subs([(q,2),(L,4)]), x)
# print roots
# print 'Solutions:'
# for z in roots:
#     ev = -a-b + a/z + b*z
#     print ev.subs([(a,2), (b,1)]).evalf()


evlist = [-a-b+a+b, -a-b+sqrt(2*a*b), -a-b,
        -a-b-sqrt(2*a*b),
        (-3*a-3*b+sqrt(a*a+10*a*b+b*b))/2,
        (-3*a-3*b-sqrt(a*a+10*a*b+b*b))/2]

print 'benchmark:'
for val in evlist:
    print val.subs([(a,2), (b,1)]).evalf()

