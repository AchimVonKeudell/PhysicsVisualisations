# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 07:41:41 2026

@author: Achim
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter

# Domain
N = 128              # number of spatial points
L = 32*np.pi         # domain length
dx = L / N
x = np.linspace(0, L, N, endpoint=False)

# Time
dt = 0.001
steps = 400000
plot_every = 500
framesnumber = int(steps/plot_every)

# Wavenumbers
k = 2*np.pi*np.fft.fftfreq(N, d=dx)
k2 = k**2
k4 = k**4

# Initial condition
u = np.cos(x/16)*(1 + np.sin(x/16))

def rhs(u):
    """Right-hand side of Kuramoto-Sivashinsky equation."""
    u_hat = np.fft.fft(u)
    ux = np.real(np.fft.ifft(1j*k*u_hat))
    uxx = np.real(np.fft.ifft(-k2*u_hat))
    uxxxx = np.real(np.fft.ifft(k4*u_hat))
    return -uxx - uxxxx - u*ux

# Storage for visualization
frames1 = []

# Time integration (RK4)
for n in range(steps):
    k1 = rhs(u)
    k2_ = rhs(u + 0.5*dt*k1)
    k3 = rhs(u + 0.5*dt*k2_)
    k4_ = rhs(u + dt*k3)
    u = u + dt*(k1 + 2*k2_ + 2*k3 + k4_) / 6

    if n % plot_every == 0:
        frames1.append(u.copy())

# Animation
fig, ax = plt.subplots()
line, = ax.plot(x, frames1[0])
ax.set_ylim(-5,5)
ax.set_xlim(0,100)
ax.set(xlabel='x',ylabel='u')
infobox = '$\partial_tu + \partial_{x}^2u + \partial_{x}^4u + 0.5(\partial_{x}u)^2 = 0$'
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 

ax.text(0.05,0.95,infobox, bbox=props,fontsize=12,verticalalignment='top',transform=ax.transAxes)
fig.suptitle('Kuramoto–Sivashinsky equation')

def update(k):
    line.set_ydata(frames1[k])
    return line,

#writer = FFMpegWriter(fps=30)
ani = FuncAnimation(fig, update, frames=framesnumber)
ani.save('kuramoto1d.gif',fps=25,dpi=300) 
#ani.save('kuramoto.mp4',writer=writer) 

