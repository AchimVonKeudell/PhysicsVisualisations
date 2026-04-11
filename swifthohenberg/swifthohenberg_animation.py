#
# Swift Hohenberg PDE
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2026
#

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Grid parameters
N = 1024
L = 400.0
dx = L / N
frameNB = 2000

# Time parameters
dt = 0.1
steps_per_frame = 5

# Control parameter
r = 0.7

# Create grid
x = np.linspace(-L/2, L/2, N)
y = np.linspace(-L/2, L/2, N)
X, Y = np.meshgrid(x, y)

# Initial condition (small noise)
u = 0.1 * np.random.randn(N, N)

# Fourier space setup
kx = 2 * np.pi * np.fft.fftfreq(N, d=dx)
ky = 2 * np.pi * np.fft.fftfreq(N, d=dx)
KX, KY = np.meshgrid(kx, ky)
k2 = KX**2 + KY**2


# Setup plot
fig, ax = plt.subplots(figsize=(6,6))

im = ax.imshow(u, cmap='plasma', origin='lower', extent=[-L/2, L/2, -L/2, L/2])

infobox = '$\partial_t u = r u - (1 + \partial_{xy}^2)^2 u -u^3$'
#infobox += 'F: '+"{:.3f}".format(F) + ', k: '
props = dict(boxstyle='round', facecolor='lightblue', alpha=1) 
ax.set(xlabel='x',ylabel='y')
ax.text(0.05,0.95,infobox, bbox=props,fontsize=10,verticalalignment='top',transform=ax.transAxes)
ax.axis('off')

def step(u,r):
    u_hat = np.fft.fft2(u)

    # Nonlinear term in real space
    nonlinear = -u**3
    nonlinear_hat = np.fft.fft2(nonlinear)
    
    # Linear operator in Fourier space
    L_hat = r - (1 - k2)**2

    # Semi-implicit time stepping
    u_hat_new = (u_hat + dt * nonlinear_hat) / (1 - dt * L_hat)

    return np.real(np.fft.ifft2(u_hat_new))

def update(frame):
    global u
    r = 0.5+frame/frameNB
    for _ in range(steps_per_frame):
        u = step(u,r)

    im.set_array(u)
    fig.suptitle('Swift–Hohenberg Pattern, r = '+ "{:.3f}".format(r))
    return [im]

anim = FuncAnimation(fig, update, frames=frameNB, interval=50)
anim.save('swifthohenberg.gif',fps=50,dpi=300) 

plt.show()