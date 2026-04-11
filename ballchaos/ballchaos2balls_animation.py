#
# Bouncing Ball inside a Circle
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2026
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection

 
filename = 'ballchaos2balls_kinetic.gif' 
framesNB = 3000

gravity = np.array([0, -0.3])  # gravity vector
radius_circle = 5
radius_ball = 0.2
restitution = 1  # bounce energy retention
dt = 0.3
x0 = np.array([2.0,1.0])
v0 = np.array([0.0,1.0])
x1 = np.array([2.0,1.0+1e-3])
v1 = np.array([0.0,1.0-1e-3])
xbase = -radius_circle


E0 = (x0[1]-xbase)*np.abs(gravity[1]) + 0.5 *np.linalg.norm(v0)**2
E1 = (x1[1]-xbase)*np.abs(gravity[1]) + 0.5 *np.linalg.norm(v1)**2


colour1 = np.array([255, 165, 0])/255
colour2 = np.array([0, 165, 255])/255

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
pos0 = x0 # np.array([0.0, 2.0])
vel0 = v0 # np.array([2.0, 0.0])
pos1 = x1 # np.array([0.0, 2.0])
vel1 = v1 # np.array([2.0, 0.0])


# Setup plot
fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.set_aspect('equal')
ax.set_axis_off()

# Draw circle boundary
circle = plt.Circle((0, 0), radius_circle, fill=False)
ax.add_patch(circle)

# Draw ball 0
ball0 = plt.Circle(pos0, radius_ball, color='red')
ax.add_patch(ball0)
trace0 = ax.plot(0,0, linestyle='dashed',lw=0.5,color='orange')[0]
trace0x = []
trace0y = []
x = [0]
y = [0]
lc0 = fade_line(x, y, colour1, 1)
trace0lc = ax.add_collection(lc0)

# Draw ball 1
ball1 = plt.Circle(pos1, radius_ball, color='blue')
ax.add_patch(ball1)
trace1 = ax.plot(0,0, linestyle='dashed',lw=0.5,color='lightblue')[0]
trace1x = []
trace1y = []
x = [0]
y = [0]
lc1 = fade_line(x, y, colour2, 1)
trace1lc = ax.add_collection(lc1)


def update(frame):
    global pos0, vel0, pos1, vel1, trace0lc, trace1lc

    # Apply gravity
    vel0 += gravity*dt

    # Update position
    pos0 += vel0*dt
        
    # Check energy conservation
    vel0th = np.sqrt(2*(E0 - (pos0[1]-xbase)*np.abs(gravity[1])))
    factor = vel0th/np.linalg.norm(vel0) 
    vel0 *= factor    
    
    # Check collision with circle boundary
    dist = np.linalg.norm(pos0)
    if dist + radius_ball >= radius_circle:
        # Normal vector
        normal = pos0 / dist

        # Reflect velocity
        vel0 = vel0 - 2 * np.dot(vel0, normal) * normal

       
        # Push ball back inside
        pos0 = normal * (radius_circle - radius_ball)

    # Update drawing
    ball0.center = pos0
    
    trace0x.append(ball0.center[0])
    trace0y.append(ball0.center[1])
    trace0lc.remove()
    x = trace0x[-1000:]
    y = trace0y[-1000:]  
    lc0 = fade_line(x, y, colour1, 1)
    trace0lc = ax.add_collection(lc0)
    
    #######################################
    
    # Apply gravity
    vel1 += gravity*dt

    # Update position
    pos1 += vel1*dt
        
    # Check energy conservation
    vel1th = np.sqrt(2*(E1 - (pos1[1]-xbase)*np.abs(gravity[1])))
    factor = vel1th/np.linalg.norm(vel1) 
    vel1 *= factor    
    
    # Check collision with circle boundary
    dist = np.linalg.norm(pos1)
    if dist + radius_ball >= radius_circle:
        # Normal vector
        normal = pos1 / dist

        # Reflect velocity
        vel1 = vel1 - 2 * np.dot(vel1, normal) * normal

       
        # Push ball back inside
        pos1 = normal * (radius_circle - radius_ball)

    # Update drawing
    ball1.center = pos1
    
    trace1x.append(ball1.center[0])
    trace1y.append(ball1.center[1])
    trace1lc.remove()
    x = trace1x[-1000:]
    y = trace1y[-1000:]  
    lc1 = fade_line(x, y, colour2, 1)
    trace1lc = ax.add_collection(lc1)
   
    fig.suptitle('2 Bouncing Balls ')#+ "{:.4f}".format(E1))
    
    # Time: ' +' s'
    
    return ball0,

ani = FuncAnimation(fig, update, frames=framesNB, interval=30, blit=True)
ani.save(filename,fps=50,dpi=300) 

plt.show()
