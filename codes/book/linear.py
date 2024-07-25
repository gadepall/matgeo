#Code by GVV Sharma
#July 22, 2024
#released under GNU GPL
#Line 


import sys                                          #for path to external scripts
sys.path.insert(0, '/sdcard/github/matgeo/codes/CoordGeo')        #path to my scripts
import numpy as np
import mpmath as mp
import numpy.linalg as LA
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

#local imports
from line.funcs import *
from triangle.funcs import *
from conics.funcs import circ_gen


#if using termux
import subprocess
import shlex
#end if


#Given points
A = np.array(([1, 2])).reshape(-1,1) 
B = 13/5*e1
C = np.array(([5, 3])).reshape(-1,1) 
D = np.array(([5, -3])).reshape(-1,1) 

#Generating Lines
x_AD = line_gen(A,D)
x_CB = line_gen(C,B)
x_CD = line_gen(C,D)

#Plotting all lines
plt.plot(x_AD[0,:],x_AD[1,:],label='$distance(AD)$')
plt.plot(x_CB[0,:],x_CB[1,:],label='$distance(CB)$')
plt.plot(x_CD[0,:],x_CD[1,:],label='$distance(CD)$')

colors = np.arange(1,5)
#Labeling the coordinates
tri_coords = np.block([A,B,C,D])
plt.scatter(tri_coords[0,:], tri_coords[1,:], c=colors)
vert_labels = ['P','A','Q','R']
for i, txt in enumerate(vert_labels):
    #plt.annotate(txt, # this is the text
    plt.annotate(f'{txt}\n({tri_coords[0,i]:.1f}, {tri_coords[1,i]:.1f})',
                 (tri_coords[0,i], tri_coords[1,i]), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(30,0), # distance from text to points (x,y)
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
plt.legend(loc='best')
'''
plt.grid() # minor
plt.axis('equal')

#if using termux
plt.savefig('chapters/11/10/4/22/figs/fig.pdf')
subprocess.run(shlex.split("termux-open chapters/11/10/4/22/figs/fig.pdf"))
#else
#plt.show()
