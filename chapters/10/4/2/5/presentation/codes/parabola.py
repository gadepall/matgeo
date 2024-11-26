import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

import sys
sys.path.insert(0, 'CoordGeo')

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import *

#setting up plot
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
len = 100
y = np.linspace(-10,10,len)

#conic parameters
V = np.array(([1,0],[0,0]))
u = -7/2*e1 - 1/2*e2
f = -60

n,c,F,O,lam,P,e = conic_param(V,u,f)
flen = parab_param(lam,P,u)
x = parab_gen(y,flen)
xStandard = np.block([[x],[y]])

#Affine conic generation
Of = O.flatten()
xActual = P@xStandard + Of[:,np.newaxis]

n = np.array([0, 1]).reshape(-1, 1)
c = 0
m,h = param_norm(n, c)

q = chord(V,u,f,m,h)
A = q[:, 0]
B = q[:, 1]

xAxis = line_norm(n, c, -50, 50)
plt.plot(xAxis[0,:], xAxis[1,:])

#plotting
plt.plot(xActual[0,:],xActual[1,:],label='$y = x^2 - 7x - 60$')

colors = np.arange(1,3)
#Labeling the coordinates
tri_coords = q
plt.scatter(tri_coords[0,:], tri_coords[1,:], c=colors)
vert_labels = ['$\\mathbf{A}$','$\\mathbf{B}$']
for i, txt in enumerate(vert_labels):
#    plt.annotate(txt, # this is the text
    plt.annotate(f'{txt}\n({tri_coords[0,i]:.0f}, {tri_coords[1,i]:.0f})',
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-20,5), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

plt.legend()
plt.grid() # minor

plt.savefig('../figs/parabola.pdf')