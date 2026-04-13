# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 15:20:58 2026

@author: Achim
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# -----------------------------
# Aizawa system parameters
# -----------------------------
a = 0.95
b = 0.7
c = 0.6
d = 3.5
e = 0.25
f = 0.1

dt = 0.01
steps = 50000  # total simulation steps
microsteps = 30

# -----------------------------
# Initialize arrays
# -----------------------------
x = np.zeros(steps)
y = np.zeros(steps)
z = np.zeros(steps)

# initial condition
x[0], y[0], z[0] = 0.1, 0.0, 0.0

# -----------------------------
# Integrate (Euler method)
# -----------------------------
for i in range(steps - 1):
    dx = (z[i] - b) * x[i] - d * y[i]
    dy = d * x[i] + (z[i] - b) * y[i]
    dz = c + a * z[i] - (z[i]**3)/3 - (x[i]**2 + y[i]**2)*(1 + e*z[i]) + f*z[i]*(x[i]**3)

    x[i+1] = x[i] + dx * dt
    y[i+1] = y[i] + dy * dt
    z[i+1] = z[i] + dz * dt

# -----------------------------
# Setup 3D figure
# -----------------------------
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111, projection='3d')
ax.grid(False)

#ax.set_facecolor("black")
#fig.patch.set_facecolor("black")

line, = ax.plot([], [], [], lw=0.6, color='blue')
line1, = ax.plot([], [], [], marker ='o',color='red',lw=0)


ax.set_xlim(np.min(x), np.max(x))
ax.set_ylim(np.min(y), np.max(y))
ax.set_zlim(np.min(z), np.max(z))

#ax.set_xlim(-30,30)
#ax.set_ylim(-30,30)
#ax.set_zlim(0,60)
ax.set(xlabel='x',ylabel='y',zlabel='z')

ax.set_title("Aizawa Attractor")

infobox = 'dx/dt = (z-b)x - dy \n'
infobox += 'dy/dt = dx + (z-b)y \n'
infobox += 'dz/dt = c + az -z$^3$/3 + \n'
infobox += '        (x$^2$+y$^2$)(1+ez) + fzx$^3$\n'
infobox += 'a = 0.95,  b = 0.7, c = 0.6,\n'
infobox += 'd=3.5, e=0.25, f=0.1'

props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 

ax.text(0,0,21,infobox, bbox=props,fontsize=8,verticalalignment='top',transform=ax.transAxes)
# -----------------------------
# Animation function
# -----------------------------
def update(frame):
    idx = frame * microsteps  # speed control
    line.set_data(x[:idx], y[:idx])
    line.set_3d_properties(z[:idx])
    line1.set_data([x[idx]], [y[idx]])
    line1.set_3d_properties([z[idx]])
    return line,

ani = FuncAnimation(
    fig,
    update,
    frames=range(1, steps // microsteps),
    interval=20,
    blit=False
)
ani.save('aizawa.gif',fps=25,dpi=300) 

plt.show()