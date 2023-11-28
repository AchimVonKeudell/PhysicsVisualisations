# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 08:26:45 2023

@author: Achim
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import numpy as np



dt = 1e-5 # time step
steps = 1500 # simulation steps
ministeps = 1000
trotation = steps*ministeps*dt # time for totation of the image around z-axis
elev = 30
azim = 15

# 'vibrational assisted'
vx = -0.002 # initial velcoity
vy = 1.8 # initial z velocity

# 'cold molecule'
v2x = -1.8 # initial velcoity
v2y = 0.002 # initial z velocity

gamma = 0.15

# Morse potentials
E1 = 3
a1 = 1
xgw = 3

E2 = 2
a2 = 2
ygw = 2

# Gaussian potentials
xtrap = xgw
ytrap = 4
trapdepth = 4
a3 = 1

def pot(x1,y1):
   potential = (E1*(1-np.exp(-a1*(x1-xgw)))**2 + E2*(1-np.exp(-a2*(y1-ygw)))**2
               - trapdepth*np.exp(-((x1-xtrap)**2+(y1-ytrap)**2)/a3))   
   return potential  

def gradpotx(x1,y1):
  # potential = -(pot(x1+dx,y1)-pot(x1,y1))/dx*f      
   potential = (-1)*(E1*2*a1*(1-np.exp(-a1*(x1-xgw)))
                            *np.exp(-a1*(x1-xgw)) 
                     +1/a3*trapdepth*(x1-xtrap)*2*
                      np.exp(-((x1-xtrap)**2+(y1-ytrap)**2)/a3))
   return potential  

def gradpoty(x1,y1):
  # potential = -(pot(x1,y1+dx)-pot(x1,y1))/dx*f      
   potential = (-1)*(E2*2*a2*(1-np.exp(-a2*(y1-ygw)))
                            *np.exp(-a2*(y1-ygw))
                     +1/a3*trapdepth*(y1-ytrap)*2*
                      np.exp(-((x1-xtrap)**2+(y1-ytrap)**2)/a3))
   return potential  

zscalemin = -2
zscalemax = 15

x = np.linspace(2.1, 6, 100)
y = np.linspace(1.5, 6, 100)

x, y = np.meshgrid(x, y)
z = pot(x,y)


# --------------------------------------------------------------------------
# Timing Simulation
# --------------------------------------------------------------------------
t = 0


x0 = 6
y0 = 2 
z0 = pot(x0,y0)

xp = x0
yp = y0
zp = z0  

x2p = x0
y2p = y0
z2p = z0  

# reserve memory for trajectorys
xt = []
yt = []
zt = []

x2t = []
y2t = []
z2t = []


# Create a 3D plot
fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(x, y, z, vmin= zscalemin, vmax = zscalemax,color='gray',shade=False,edgecolor='None',alpha=0.5)
contour_levels = np.linspace(-2,10,30)
for h in contour_levels:
  cset = ax.contour(x, y, z, [h], vmin = 1, vmax=zscalemax, offset=h,colors='green',linewidths=0.5,linestyles='solid')


ax.view_init(elev, azim)
ax.set_zlim(zscalemin,zscalemax)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_xlabel('x distance to surface')
ax.set_ylabel('y distance in molecule')
ax.set_zlabel('potential energy')


# trajectory line
line1 = ax.plot3D(xt,yt,zt,color='r',lw=1,linestyle='dashed')[0]
line1p = ax.plot3D(xp,yp,zp,color='r',marker='o',markersize=2,lw=0,label='hot molecule')[0]
line2 = ax.plot3D(x2t,y2t,z2t,color='b',lw=1,linestyle='dashed')[0]
line2p = ax.plot3D(x2p,y2p,z2p,color='b',marker='o',markersize=2,lw=0,label='cold molecule')[0]

ax.legend(fontsize=10)
# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
    global t,dt
    global xt,yt,zt
    global xp,yp,zp
    global vx,vy
    global x2t,y2t,z2t
    global x2p,y2p,z2p
    global v2x,v2y
    
    for i in range(ministeps):
      t += dt
      vx += gradpotx(xp,yp) * dt - gamma * vx * dt     
      vy += gradpoty(xp,yp) * dt - gamma * vy * dt     
      xp += dt * vx 
      yp += dt * vy 
      zp = pot(xp,yp)
   
      v2x += gradpotx(x2p,y2p) * dt - gamma * v2x * dt
      v2y += gradpoty(x2p,y2p) * dt - gamma * v2y * dt     
      x2p += dt * v2x      
      y2p += dt * v2y 
      z2p = pot(x2p,y2p)
   
    # Update the new positions in vector
    xt.append(xp)
    yt.append(yp)
    zt.append(zp)
    line1.set_xdata(xt)
    line1.set_ydata(yt) 
    line1.set_3d_properties(zt) 
   
    xpn = np.array([xp])
    ypn = np.array([yp])
    zpn = np.array([zp])
    line1p.set_xdata(xpn)
    line1p.set_ydata(ypn) 
    line1p.set_3d_properties(zpn) 

    x2t.append(x2p)
    y2t.append(y2p)
    z2t.append(z2p)
    line2.set_xdata(x2t)
    line2.set_ydata(y2t) 
    line2.set_3d_properties(z2t) 
   
    x2pn = np.array([x2p])
    y2pn = np.array([y2p])
    z2pn = np.array([z2p])
    line2p.set_xdata(x2pn)
    line2p.set_ydata(y2pn) 
    line2p.set_3d_properties(z2pn) 

    ax.view_init(elev,azim+t/trotation*360)  
    
# ----------------------------------------
# Create the complete time sequence
# teh total time span is tend/dtplot
# returns the rate in
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save('sticking_hotcold.mp4',fps=25,dpi=300)
