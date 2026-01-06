#
# Phytagoras Tree
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2026
#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.animation  as animation





mode = 'TRLcorner'
mode = 'TRLcorner'

if mode == 'TRLcorner':
    name = 'Phytagoras tree, Top corners'
    animationname = 'phytagoras_TRLcorners.gif'
    cornermode = ''
    levels =14       # number of generations
    angledeg = 30               # branching angle
    base_length = 0.2            # base square size
    sweep = 20
    xmin=-2.2
    xmax=1.3
    ymin=-0.5
    ymax=3

if mode == 'childcorner':
    name = 'Phytagoras tree, , Child corners'
    animationname = 'phytagoras_childcorner.gif'
    cornermode = 'childcorner'
    levels =14       # number of generations
    angledeg = 30               # branching angle
    base_length = 0.2            # base square size
    sweep = 20
    xmin=-2.2
    xmax=1.3
    ymin=-0.5
    ymax=3



# ------------------ PARAMETERS ------------------

angle = np.radians(angledeg)           
scale_left = np.cos(angle)  # scaling factors
scale_right = np.sin(angle)


def square_vertices(x, y, L, theta):
    """
    Return the 4 vertices of a square given bottom-left corner,
    side length L, and rotation angle theta.
    """
    R = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta),  np.cos(theta)]
    ])

    corners = np.array([
        [0, 0],
        [L, 0],
        [L, 3*L],
        [0, 3*L]
    ])

    rotated = corners @ R.T
    return rotated + np.array([x, y])


def draw_square(ax, verts, color,alpha):
    square = Polygon(verts, closed=True, facecolor='white', edgecolor='white',alpha=alpha)
    ax.add_patch(square)
    square = Polygon(verts, closed=True, facecolor=color, edgecolor=color,alpha=alpha)
    ax.add_patch(square)


# ------------------ STORAGE ------------------
# Each square: (x, y, L, theta)
generations = []
generations.append([(0.0, 0.0, base_length, 0.0)])

# ------------------ GENERATE TREE ------------------
for gen in range(levels - 1):
    next_gen = []
    for x, y, L, theta in generations[gen]:

        # Parent square vertices
        verts = square_vertices(x, y, L, theta)

        # Top-left and top-right corners
        TL = verts[3]
        TR = verts[2]

        # Left child
        L_left = L * scale_left
        theta_left = theta + angle
        next_gen.append((TL[0], TL[1], L_left, theta_left))
        
        # Right child
        L_right = L * scale_right
        theta_right = theta - (np.pi / 2 - angle)
        if cornermode == 'childcorner':
          vertsleftchild = square_vertices(TL[0], TL[1], L_left, theta_left)
          BR = vertsleftchild[1]
          next_gen.append((BR[0], BR[1], L_right, theta_right))
        else:
          next_gen.append((TR[0], TR[1], L_right, theta_right))  

    generations.append(next_gen)


fig, ax = plt.subplots(figsize=(5, 5))

ax.set(xlim=(xmin,xmax))
ax.set(ylim=(ymin,ymax))
ax.set_aspect("equal")

ax.axis("off")
# info box
infobox = name+'\n'
infobox +='generations: ' + "{:.0f}".format(levels)+'\n'
infobox +='angle: ' + "{:.0f}".format(angledeg)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)


# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
    
  kindex = int(k/sweep)  
  print(kindex,' von ',levels,' branches ',len(generations[kindex])) 
    
  gen = generations[kindex]    
  color = plt.cm.plasma(kindex / levels)
  for x, y, L, theta in gen:
        verts = square_vertices(x, y, L, theta)
        draw_square(ax, verts, color,alpha=(k-kindex*sweep)/(sweep))

  fig.suptitle('Phytagoras Tree, Generation: '+"{:.0f}".format(kindex))


anim = animation.FuncAnimation(fig,animate,interval=1,frames=levels*sweep)
anim.save(animationname,fps=25,dpi=300)       
    

