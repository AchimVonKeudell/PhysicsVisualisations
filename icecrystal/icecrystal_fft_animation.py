# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 16:27:23 2024

@author: Achim
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation


n = 3000 # grid size
freezingthreshold = 0.5
animationfile = 'icegrowth.mp4'
iterations = 2000

# Preferred direction on a hexagonal grid
def getSmoother(n):
    smoother = np.zeros([n, n])
    smoother[0, 0] = 1
    smoother[0, 1] = 1
    smoother[1, 0] = 1
    smoother[-1, 0] = 1
    smoother[0, -1] = 1
    smoother[1, -1] = 1
    smoother[-1, 1] = 1
    return smoother

# Gaussian Kernel to account for diffusion
def getG(n):
    G = np.zeros([n, n])
    for i in range(n):
        for j in range(n):
            i2 = i
            j2 = j
            if i2 > n/2:
                i2 -= n
            if j2 > n/2:
                j2 -= n
            # we're on a hexagonal grid, so this
            # is the x displacement for two tiles
            dx = i2 + 0.5 * j2
            # and the y displacement
            dy = j2 * math.sqrt(3)/2
            r = math.sqrt(dx*dx+dy*dy)
            G[i, j] = math.exp(-r*r/8)
    G = G / np.sum(np.sum(G))
    return G

# Design Plot 
fig, ax = plt.subplots(1,1,figsize=(5,5))
A = np.zeros([n, n])
im = ax.imshow(A,cmap='Greys',vmin=0,vmax=1,label='ice')
imrho = ax.imshow(A,cmap='Blues',vmin=0,vmax=1,alpha=0.5,label='water')
#ax.set(xlabel="x",ylabel="h")
#ax.legend(loc=1,fontsize=6)
ax.set_xticks([])
ax.set_yticks([])
fig.colorbar(imrho,ax=ax,label='water density')

# setup info box
infobox = ''
infobox += 'freezing threshold: ' + "{:.2f}".format(freezingthreshold) + '' 
#infobox += 'facet angle spread: ' + "{:.0f}".format(anglewidth/np.pi*180)   
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.98,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)

# Initialize Simulation
A = np.zeros([n, n])
A[n//2, n//2] = 1
G = getG(n)
smoother = getSmoother(n)
rho = np.ones([n, n])

sf = np.fft.fft2(smoother)
gf = np.fft.fft2(G)
  

def animate(k):
  
  global A,rho,n  
  print(k)
  
  # microsteps
  for i in range(5):  
     B = A.copy()
     B[B < 1] = 0
     B[B > 1] = 1

     # propagate freezing by convoluting with crystal pattern
     A2 = np.fft.ifft2(np.multiply(sf, np.fft.fft2(B))).real
     # decide freezing with threshold
     A2[A2 > freezingthreshold] = 1
     A2[A2 < freezingthreshold] = 0
     
     # add average water
     A += A2 * (rho * 0.5)
     
     # upper boundary of water
     A[A > 5] = 5
     
     # remove water close to ice crystal
     rho[A > 0.5] = 0
     
     # diffuse water by convoluting with Gaussian Kernel
     rho = np.fft.ifft2(np.multiply(gf, np.fft.fft2(rho))).real
     
     # normalize rho and add randm pertrubation  
     rho = rho / np.mean(np.mean(rho)) + np.random.randn(n, n) * 0.02
  
  B = A * 0.2
  n = len(A)
  B[B < 0.01] = 0

  # This handles the hexagonal grid, somewhat approximately:
  B = np.array(
        [B[i, range((n // 2 - i // 2), (n - i // 2), 1)] for i in range(n // 4, (3 * n) // 4, 1) if (i % 5) != 0])
  rhomat = np.array(
        [rho[i, range((n // 2 - i // 2), (n - i // 2), 1)] for i in range(n // 4, (3 * n) // 4, 1) if (i % 5) != 0])

  #def red(x):
  #      return np.uint(math.exp(-x * 4.0) * 255)

  #def green(x):
  #      return np.uint(math.exp(-x) * 255)

  #def blue(x):
  #      return np.uint(math.exp(-x / 4.0) * 255)

  #color = [[[red(x), green(x), blue(x), 255] for x in y] for y in B]
  #im.set_array(color)
  im.set_array(B)
  imrho.set_array(rhomat)


anim = animation.FuncAnimation(fig,animate,interval=1,frames=iterations)
anim.save(animationfile,fps=25,dpi=300)