# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 12:36:39 2026

@author: Achim
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from numba import njit

model = 'chaos'
model = 'pattern'
model = 'transition'

# Grid parameters
Nx, Ny = 100, 100
dx = dy = 1.0
dt = 0.01

# Initialize field with small random noise
u = 0.1 * np.random.randn(Nx, Ny)

@njit
def laplacian(u, dx, dy):
    Nx, Ny = u.shape
    lap = np.zeros_like(u)
    for i in range(Nx):
        for j in range(Ny):
            lap[i, j] = (
                u[(i+1)%Nx, j] + u[(i-1)%Nx, j]
                + u[i, (j+1)%Ny] + u[i, (j-1)%Ny]
                - 4*u[i, j]
            ) / dx**2
    return lap

@njit
def biharmonic(u, dx, dy):
    return laplacian(laplacian(u, dx, dy), dx, dy)

@njit
def grad_sq(u, dx, dy):
    Nx, Ny = u.shape
    g2 = np.zeros_like(u)
    for i in range(Nx):
        for j in range(Ny):
            ux = (u[(i+1)%Nx, j] - u[(i-1)%Nx, j]) / (2*dx)
            uy = (u[i, (j+1)%Ny] - u[i, (j-1)%Ny]) / (2*dy)
            g2[i, j] = ux**2 + uy**2
    return g2

@njit
def step(u, dt, dx, dy, a):
    lap = laplacian(u, dx, dy)
    bih = biharmonic(u, dx, dy)
    g2 = grad_sq(u, dx, dy)
    
    #return u + dt * (-lap - bih - 0.5 * g2)
    return u + dt * (-a*u -lap - bih + 0.5 * g2)

# Set up plot
fig, ax = plt.subplots()
im = ax.imshow(u, cmap='inferno', animated=True, vmin=-10,vmax=10)
#infobox = '$\partial_tu + \partial_{xy}^2 u + \partial_{xy}^4 u + 0.5(\partial_{xy} u)^2 = 0$'
infobox = '$\partial_tu + a u + \partial_{xy}^2 u + \partial_{xy}^4 u - 0.5(\partial_{xy} u)^2 = 0$'
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8) 
ax.set(xlabel='x',ylabel='y')
ax.text(0.05,0.95,infobox, bbox=props,fontsize=8,verticalalignment='top',transform=ax.transAxes)
plt.colorbar(im)

def update(frame):
    global u
    if model == 'chaos':
        a = 0
    if model == 'pattern':
        a = 0.225
    if model == 'transition':    
        a = 0.225*frame/2000
    for _ in range(50):  # multiple steps per frame for smoother evolution
        u = step(u, dt, dx, dy, a)
        u = u - np.mean(u)
    im.set_array(u)
    fig.suptitle('Kuramoto–Sivashinsky equation, a='+"{:.3f}".format(a))

    print(frame)
    return [im]

ani = FuncAnimation(fig, update, frames=2000, interval=50, blit=True)
ani.save('kuramoto3dc.gif',fps=25,dpi=300) 

plt.show()