# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 08:35:26 2026

@author: Achim
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection

# Parameters
gravity = np.array([0, -0.3])  # gravity vector
radius_circle = 5
radius_ball = 0.2
restitution = 1  # bounce energy retention
dt = 0.3
x0 = np.array([2.0,1.0])
v0 = np.array([0.0,1.0])
xbase = -radius_circle


E0 = (x0[1]-xbase)*np.abs(gravity[1]) + 0.5 *np.linalg.norm(v0)**2
print(E0)

colour1 = np.array([255, 165, 0])/255
colour2 = np.array([255, 165, 0])/255

def fade_line(x, y, colour, alpha):
    Npts = len(x)
    alpha = 1

    # create colours array
    colours = np.zeros((Npts, 4))
    colours[:, 0:3] = colour
    colours[:, 3] = alpha*np.linspace(0, 1, Npts)

    # N-1 segments for N points
    # (x, y) start point and (x2, y2) end point
    # in 3D array each segment 2D inside
    segments = np.zeros((Npts-1, 2, 2))
    # segment start values - slice of last value
    segments[:, 0, 0] = x[:-1]
    segments[:, 0, 1] = y[:-1]
    # segements end values - slice of first value
    segments[:, 1, 0] = x[1:]
    segments[:, 1, 1] = y[1:]
    
    lc = LineCollection(segments, color=colours)
    return lc



# Initial state
pos = x0 # np.array([0.0, 2.0])
vel = v0 # np.array([2.0, 0.0])

# Setup plot
fig, ax = plt.subplots()
ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.set_aspect('equal')
ax.set_axis_off()

# Draw circle boundary
circle = plt.Circle((0, 0), radius_circle, fill=False)
ax.add_patch(circle)

# Draw ball
ball = plt.Circle(pos, radius_ball, color='red')
ax.add_patch(ball)

trace = ax.plot(0,0, linestyle='dashed',lw=0.5,color='orange')[0]

tracex = []
tracey = []
x = [0]
y = [0]
lc = fade_line(x, y, colour1, 1)
tracelc = ax.add_collection(lc)
usefadeline = True


def update(frame):
    global pos, vel, tracelc

    # Apply gravity
    vel += gravity*dt

    # Update position
    pos += vel*dt
    
    # Check energy conservation
    #Ekin = 0.5 *np.linalg.norm(vel)**2
    #Ekinthe = (x0[1]-pos[1])*np.abs(gravity[1])
    #print(Ekin,Ekinthe,x0[1],pos[1],np.abs(gravity[1]))
    #factor = Ekin/Ekinthe
    #vel *= np.sqrt(factor)
    #print(factor)
    E1 = (pos[1]-xbase)*np.abs(gravity[1])+0.5*np.sum(np.square(vel))
    factor=E0/E1
    vel *= np.sqrt(factor)
    
    # Check collision with circle boundary
    dist = np.linalg.norm(pos)
    if dist + radius_ball >= radius_circle:
        # Normal vector
        normal = pos / dist

        # Reflect velocity
        vel = vel - 2 * np.dot(vel, normal) * normal

        # Apply restitution
        vel *= restitution

        # Push ball back inside
        pos = normal * (radius_circle - radius_ball)

    # Update drawing
    ball.center = pos
    
    tracex.append(ball.center[0])
    tracey.append(ball.center[1])
    if usefadeline:
      # create fade line
      tracelc.remove()
      x = tracex[-1000:]
      y = tracey[-1000:]  
      lc = fade_line(x, y, colour1, 1)
      tracelc = ax.add_collection(lc)
    else: 
      trace.set_xdata(tracex)
      trace.set_ydata(tracey)
 
    fig.suptitle('Bouncing Ball ')#+ "{:.4f}".format(E1))
    
    # Time: ' +' s'
    
    return ball,

ani = FuncAnimation(fig, update, frames=4000, interval=30, blit=True)
ani.save('ballchaos.gif',fps=50,dpi=300) 


plt.title("Ball Bouncing Inside a Circle")
plt.show()
