#
# Theodorus spiral
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2026
#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Polygon
import matplotlib.cm as cm

# ----------------------------
# Build Theodorus spiral points
# ----------------------------
def generate_theodorus(n_steps):
    points = [np.array([0.0, 0.0]), np.array([1.0, 0.0])]
    angle = 0.0

    for i in range(2, n_steps + 1):
        radius = np.sqrt(i)
        angle += np.arctan(1 / np.sqrt(i - 1))
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        points.append(np.array([x, y]))

    return points


# ----------------------------
# Create triangles from points
# ----------------------------
def build_triangles(points):
    triangles = []
    origin = points[0]

    for i in range(1, len(points)):
        triangles.append([origin, points[i - 1], points[i]])

    return triangles


# ----------------------------
# Parameters
# ----------------------------
N = 80  # number of triangles
fadeNB = 15
points = generate_theodorus(N)
triangles = build_triangles(points)

# Colormap
cmap = cm.get_cmap("viridis", N+1)

# ----------------------------
# Plot setup
# ----------------------------
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_aspect("equal")
ax.set_xlim(-1-np.sqrt(N), np.sqrt(N) + 1)
ax.set_ylim(-1-np.sqrt(N), np.sqrt(N) + 1)
ax.axis('off')

patches = []


# ----------------------------
# Animation function
# ----------------------------
def update(frame):
    triangleIDX = int(frame/fadeNB)
    alphaID = (frame-triangleIDX*fadeNB)/fadeNB
    if triangleIDX < len(triangles):
        tri = triangles[triangleIDX]
        polygon = Polygon(
            tri,
            closed=True,
            facecolor=cmap(triangleIDX),
            edgecolor="black",
            alpha=alphaID,
        )
        ax.add_patch(polygon)
        if len(patches)>1:
          if  patches[-2].get_alpha()<(1-1/fadeNB):
             patches[-2].set_alpha(0)
        patches.append(polygon)
        
    fig.suptitle('Theodorus Spiral, length: '+"{:.0f}".format(triangleIDX)+'$^{1/2}$') #+ )
    
    return patches


# ----------------------------
# Run animation
# ----------------------------
ani = FuncAnimation(
    fig,
    update,
    frames=len(triangles)*fadeNB,
    interval=400,
    blit=False,
    repeat=False,
)

ani.save('theodorus.gif',fps=25,dpi=300) 


