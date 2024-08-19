#Program to plot an ellipse 
#Code by GVV Sharma
#August 8, 2020
#Revised July 31, 2024
#Revised August 16, 2024

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

import sys                                          #for path to external scripts
sys.path.insert(0, '/sdcard/github/matgeo/codes/CoordGeo')        #path to my scripts


#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import *

#if using termux
import subprocess
import shlex
#end if

#setting up plot
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
len = 100
y = np.linspace(-5,5,len)

#Ellipse parameters
V = np.array(([25,0],[0,4]))
u = np.array(([0,0])).reshape(-1,1)
f = -4*25
n,c,F,O,lam,P = conic_param(V,u,f)
ab = ellipse_param(V,u,f)

#Eigenvalues and eigenvectors
print(lam, P)
xStandard= ellipse_gen(ab[0],ab[1])

#Directrix
k1 = -1
k2 = 1

#Latus rectum
cl = (n.T@F).flatten()

#Generating lines
x_A = line_norm(n,c[0],k1,k2)#directrix
x_B = line_norm(n,cl[0],k1,k2)#latus rectum
x_C = line_norm(n,c[1],k1,k2)#directrix
x_D = line_norm(n,cl[1],k1,k2)#latus rectum
#Affine conic generation
Of = O.flatten()
#xActual = P@xStandard + Of[:,np.newaxis]
xActual = xStandard 

#plotting
plt.plot(xActual[0,:],xActual[1,:],label='Ellipse')
plt.plot(x_A[0,:],x_A[1,:],label='Directrix')
plt.plot(x_B[0,:],x_B[1,:],label='Latus Rectum')
plt.plot(x_C[0,:],x_C[1,:])
plt.plot(x_D[0,:],x_D[1,:])
#
colors = np.arange(1,4)
#Labeling the coordinates
tri_coords = np.block([O,F])
plt.scatter(tri_coords[0,:], tri_coords[1,:], c=colors)
vert_labels = ['$\mathbf{O}$','$\mathbf{F}_1$','$\mathbf{F}_2$']
for i, txt in enumerate(vert_labels):
#    plt.annotate(txt, # this is the text
    plt.annotate(f'{txt}\n({tri_coords[0,i]:.2f}, {tri_coords[1,i]:.2f})',
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-20,5), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center

# use set_position
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['left'].set_position('zero')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_position('zero')
'''
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
plt.xlabel('$x$')
plt.ylabel('$y$')
'''
plt.legend(loc='best')
plt.grid() # minor
plt.axis('equal')

#if using termux
plt.savefig('chapters/11/11/3/2/figs/fig.pdf')
subprocess.run(shlex.split("termux-open chapters/11/11/3/2/figs/fig.pdf"))
#else
#plt.show()