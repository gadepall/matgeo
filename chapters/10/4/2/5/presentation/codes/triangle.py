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

#vertices
A = np.array([0, 0]).reshape(-1, 1)
B = np.array([12, 0]).reshape(-1, 1)
C = np.array([0, 7]).reshape(-1, 1)

xAB = line_gen(A, B)
plt.plot(xAB[0,:], xAB[1,:], label="AB")
xBC = line_gen(B, C)
plt.plot(xBC[0,:], xBC[1,:], label="BC")
xAC = line_gen(A, C)
plt.plot(xAC[0,:], xAC[1,:], label="AC")

points = np.hstack([A, B, C])
plt.scatter(points[0,:], points[1,:])
verts = [A, B, C]
labels = ['A', 'B', 'C']
for i in range(len(verts)):
    x, y = verts[i][0, 0], verts[i][1, 0]
    plt.annotate(f"{labels[i]}\n({x}, {y})",
                 (x, y),
                 textcoords="offset points",
                 xytext=(20, -10),
                 ha="center")

plt.legend()
plt.axis('equal')
plt.grid()
plt.savefig('../figs/triangle.pdf')