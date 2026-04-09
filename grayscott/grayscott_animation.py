# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 23:17:20 2026

@author: Achim
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

model = 'flower'
model = 'mazes'
model = 'mitosis'
model = 'soliton'


if model == 'flower':
  F, k = 0.055, 0.062
  filename = 'grayscott_flower.gif'
  simtitle = 'Flower pattern'
  N = 1000
if model == 'mitosis':
  F, k = 0.028, 0.062
  filename = 'grayscott_mitosis.gif'
  simtitle = 'Mitosis pattern'
  N = 1000
if model == 'soliton':
  F, k = 0.03, 0.06
  filename = 'grayscott_soliton.gif'
  simtitle = 'Soliton pattern'
  N = 1000
if model == 'mazes':
  F, k = 0.029, 0.057
  filename = 'grayscott_mazes.gif'
  simtitle = 'Mazes pattern'
  N = 1000


# Gray-Scott parameters
Du, Dv = 0.16, 0.08


# Grid
u = np.ones((N, N))
v = np.zeros((N, N))

# Initial disturbance
r = 20
u[N//2-r:N//2+r, N//2-r:N//2+r] = 0.50
v[N//2-r:N//2+r, N//2-r:N//2+r] = 0.25

# Laplacian function
def laplacian(Z):
    return (
        -4 * Z
        + np.roll(Z, (0, -1), (0, 1))
        + np.roll(Z, (0, 1), (0, 1))
        + np.roll(Z, (-1, 0), (0, 1))
        + np.roll(Z, (1, 0), (0, 1))
    )

# Simulation parameters
dt = 1.0
total_time = 2000.0
steps = int(total_time / dt)

# Setup plot
fig, ax = plt.subplots()
im = ax.imshow(v, cmap='inferno', interpolation='bilinear')
infobox = '$\partial_tu = D_u\partial_{xy}^2 u - uv^2+ F(1-u)$ \n'
infobox += '$\partial_tv = D_v\partial_{xy}^2 v + uv^2 - (F+k)v$ \n'
infobox += 'F: '+"{:.3f}".format(F) + ', k: '+"{:.3f}".format(k)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8) 
ax.set(xlabel='x',ylabel='y')
ax.text(0.05,0.95,infobox, bbox=props,fontsize=6,verticalalignment='top',transform=ax.transAxes)


ax.set_title("Gray-Scott Model, "+simtitle)
ax.axis('off')

# Update function
def update(frame):
    global u, v
    for _ in range(10):  # sub-steps for stability
        Lu = laplacian(u)
        Lv = laplacian(v)

        uvv = u * v * v
        u += (Du * Lu - uvv + F * (1 - u)) * dt
        v += (Dv * Lv + uvv - (F + k) * v) * dt

    im.set_data(v)
    return [im]


# Animation
ani = animation.FuncAnimation(fig, update, frames=steps, interval=100)
ani.save(filename,fps=50,dpi=300) 

plt.show()