#Program to plot  the tangent of a parabola
#Code by GVV Sharma
#Released under GNU GPL
#August 10, 2020
#Revised July 31, 2024

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
x = np.linspace(-30,30,len)
D = np.array(([18,1944/625])).reshape(-1,1)
E = 18*e1

'''
#conic parameters
V = np.array(([0,0],[0,1]))
#u = -8*e1
P = rotmat(-np.pi/2)
u = -625/12*P@e1
print( -625/12*P@e2)
f = 0

n,c,F,O,lam,P,e = conic_param(V,u,f)
#print(n,c,F)
#print(lam,P)

#Eigenvalues and eigenvectors

#flen = parab_param(lam,P,u)
print(flen,e)
'''
flen = 625/6

#Standard parabola generation
y = parab_gen(x,flen)

#Cable
x_A = line_gen(D,E)


#Directrix
k1 = -8
k2 = 8

#Latus rectum
#cl = (n.T@F).flatten()

P = rotmat(np.pi/2)
#Of = O.flatten()
#F = P@F
#Generating lines
#x_A = P@line_norm(n,c,k1,k2)+ Of[:,np.newaxis]#directrix
#x_B = P@line_norm(n,cl[0],k1,k2)+ Of[:,np.newaxis]#latus rectum
#print(n,c)
xStandard =np.block([[x],[y]])

#Affine conic generation
#xActual = P@xStandard + Of[:,np.newaxis]

#plotting
#plt.plot(xActual[0,:],xActual[1,:],label='Parabola',color='r')
plt.plot(xStandard[0,:],xStandard[1,:],label='Cable',color='r')
plt.plot(x_A[0,:],x_A[1,:],label='Wire')
#plt.plot(x_B[0,:],x_B[1,:],label='Latus Rectum')
#
colors = np.arange(1,3)
#Labeling the coordinates
tri_coords = np.block([D,E])
plt.scatter(tri_coords[0,:], tri_coords[1,:], c=colors)
#vert_labels = ['$\mathbf{O}$','$\mathbf{F}$']
vert_labels = ['$\mathbf{D}$','$\mathbf{E}$']
for i, txt in enumerate(vert_labels):
    plt.annotate(txt, # this is the text
#    plt.annotate(f'{txt}\n({tri_coords[0,i]:.2f}, {tri_coords[1,i]:.2f})',
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
plt.savefig('chapters/11/11/5/3/figs/fig-temp.pdf')
subprocess.run(shlex.split("termux-open chapters/11/11/5/3/figs/fig-temp.pdf"))
#else
#plt.show()
