# plot.py
import numpy as np
import matplotlib.pyplot as plt

# Line definitions (exact)
P1 = np.array([0.0, 0.0, 1.0])    # body diagonal start
D1 = np.array([1.0, 1.0, -1.0])   # body diagonal direction -> end (1,1,0)
P2 = np.array([0.0, 0.0, 0.0])    # face diagonal start
D2 = np.array([1.0, 0.0, 1.0])    # face diagonal direction -> end (1,0,1)

# exact shortest distance using cross product
cross = np.cross(D1, D2)
distance = abs(np.dot(P2 - P1, cross)) / np.linalg.norm(cross)
print("Shortest distance (NumPy):", distance)

# compute closest points (solve 2x2 linear system)
w0 = P1 - P2
A = np.dot(D1, D1)
B = np.dot(D1, D2)
C = np.dot(D2, D2)
rhs = np.array([-np.dot(D1, w0), -np.dot(D2, w0)])
M = np.array([[A, -B],[B, -C]])
s, t = np.linalg.solve(M, rhs)
closest1 = P1 + s * D1
closest2 = P2 + t * D2
print("Closest points:", closest1, closest2)
print("Check distance:", np.linalg.norm(closest1 - closest2))

# plotting
def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()
    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])
    max_range = max(x_range, y_range, z_range)
    x_mid = np.mean(x_limits)
    y_mid = np.mean(y_limits)
    z_mid = np.mean(z_limits)
    half = max_range / 2
    ax.set_xlim(x_mid - half, x_mid + half)
    ax.set_ylim(y_mid - half, y_mid + half)
    ax.set_zlim(z_mid - half, z_mid + half)

fig = plt.figure(figsize=(7,7))
ax = fig.add_subplot(111, projection='3d')

vertices = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                     [0,0,1],[1,0,1],[1,1,1],[0,1,1]])
edges = [(0,1),(1,2),(2,3),(3,0),
         (4,5),(5,6),(6,7),(7,4),
         (0,4),(1,5),(2,6),(3,7)]
for e in edges:
    a = vertices[e[0]]
    b = vertices[e[1]]
    ax.plot([a[0],b[0]],[a[1],b[1]],[a[2],b[2]], color='gray', alpha=0.6)

# diagonals as precise lines from vertex to vertex (no arrowheads)
ax.plot([P1[0], P1[0]+D1[0]],[P1[1], P1[1]+D1[1]],[P1[2], P1[2]+D1[2]],
        color='blue', linewidth=2, label='Body diagonal (0,0,1)→(1,1,0)')
ax.plot([P2[0], P2[0]+D2[0]],[P2[1], P2[1]+D2[1]],[P2[2], P2[2]+D2[2]],
        color='red', linewidth=2, label='Face diagonal (0,0,0)→(1,0,1)')

# exact perpendicular segment
ax.plot([closest1[0], closest2[0]],
        [closest1[1], closest2[1]],
        [closest1[2], closest2[2]],
        color='green', linewidth=2, label=f'Perpendicular (dist={distance:.4f})')

# markers for vertices and closest points
ax.scatter(*P1, color='black', s=30)
ax.scatter(*(P1 + D1), color='black', s=30)
ax.scatter(*P2, color='black', s=30)
ax.scatter(*(P2 + D2), color='black', s=30)
ax.scatter(*closest1, color='blue', s=50)
ax.scatter(*closest2, color='red', s=50)

ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
ax.set_title('Cube — body & face diagonals')
ax.legend(loc='upper left')
set_axes_equal(ax)

plt.tight_layout()
plt.savefig("cube_lines.png", dpi=300, bbox_inches='tight')
plt.show()

