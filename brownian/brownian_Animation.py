#
# Brownian Motion
# Examples for physics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#

from numba import jit
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
from matplotlib.collections import LineCollection

# ------------------------------------------------------------------------
# parameter
# ------------------------------------------------------------------------ 

Natoms = 2000  
L = 1 # size of volume in m
Matom = 4E-3/6E23 # helium mass
Ratom = 0.007 # size of helium atom in arbitrary units
k = 1.38E-23 
T = 50 # temperature

dt = 5E-6
dtplot = 5e-6
simulationtime = 0.01
steps = int(simulationtime/dtplot)

NBbins = 50
xbin1 = np.zeros(NBbins)
xscale = np.linspace(0,L,NBbins)
xbin1 = np.zeros(NBbins)
xbin2 = np.zeros(NBbins)

poslist = []
plist = []
mlist = []
rlist = []

# ----------------------------------------------------------------------------
# Initialize atoms
# ----------------------------------------------------------------------------
Lmin = 1.1*Ratom
Lmax = L-Lmin

for i in range(Natoms+1):
    x = Lmin+(Lmax-Lmin)*np.random.rand()
    y = Lmin+(Lmax-Lmin)*np.random.rand()
    r = Ratom
    mass = Matom*r**3/Ratom**3
    pavg = np.sqrt(2.*mass*1*k*T) # average kinetic energy p**2/(2mass) = (2/2)kT
    
    phi = 2*np.pi*np.random.rand()
    px = pavg*np.cos(phi)
    py = pavg*np.sin(phi)

    poslist.append((x,y))
    plist.append((px,py))
    mlist.append(mass)
    rlist.append(r)
    
# additional atom
x = L/2
y = L/2
r = Ratom
mass = Matom*r**3/Ratom**3*50
pavg = 0 # np.sqrt(2.*mass*1*k*T) # average kinetic energy p**2/(2mass) = (2/2)kT
    
phi = 2*np.pi*np.random.rand()
px = pavg*np.cos(phi)
py = pavg*np.sin(phi)

poslist.append((x,y))
plist.append((px,py))
mlist.append(mass)
rlist.append(r)    
    
pos = np.array(poslist)
p = np.array(plist)
radius = np.array(rlist)

# ----------------------------------------------------------------------------
# Initialize plots
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(1,1,figsize=(6,6))
fig.tight_layout(pad=3)

scatter1 = ax.plot(pos[:Natoms,0],pos[:Natoms,1],marker='o',markersize=2,color='blue',lw=0)[0]
scatter2 = ax.plot(pos[Natoms+1,0],pos[Natoms+1,1],marker='o',markersize=8,color='red',lw=0,label='test particle')[0]
ax.set(xlabel='position',ylabel='position',xlim=(Lmin,Lmax),ylim=(Lmin,Lmax))
ax.set_aspect('equal')

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

def circle(t):
   angle = np.linspace(0,2*np.pi,360) 
   r = np.sqrt(1*t)
   x = L/2 + np.cos(angle)*r 
   y = L/2 + np.sin(angle)*r 
   return x,y

tracex = []
tracey = []
x = [0]
y = [0]
colour = np.array([255, 0, 0])/255
lc = fade_line(x, y, colour, 1)
tracelc = ax.add_collection(lc)

x,y = circle(0)
circleline = ax.plot(x,y,color='orange',linestyle='dashed',lw=1,label='statistical travel distance')[0]

ax.legend(loc=1)

t = 0.0
Nsteps = 0
time = 0

# initial half-step
pos[:,0] += p[:,0]/mlist*dt/2   
pos[:,1] += p[:,1]/mlist*dt/2   
     
# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------

@jit(nopython = True)
def findhitlist(pos):
     hitlist = []
     # find collisions
     for i in range(Natoms+2):
        for j in range(i,Natoms+2):
           d2 = np.linalg.norm(pos[j]-pos[i])**2 
           if d2 < (2*Ratom)**2: 
               if i!=j:
                   hitlist.append((i,j)) 
     return hitlist              
     
def animate(k):
  global pos,p,t,time
  global tracelc

  print('time: '+ "{:.2f}".format(t/1e-6)+' microseconds of '+"{:.0f}".format(simulationtime/1e-6)+' microseconds' )
  t += dt
  while time<=dtplot:
     time += dt   

     pos[:,0] += p[:,0]/mlist*dt   
     pos[:,1] += p[:,1]/mlist*dt   
     
     hitlist = findhitlist(pos)
     
     hits = np.array(hitlist)           
  
     # change momentum
     for i in range(len(hits)):
         ii = hits[i,0]
         jj = hits[i,1]         
    
         m1 = mlist[ii]
         m2 = mlist[jj]
         
         mtot = m1+m2    
    
         r1 = pos[ii]
         r2 = pos[jj]
         v1 = p[ii]/m1
         v2 = p[jj]/m2
         d = np.linalg.norm(r1 - r2)**2
         u1 = v1 - 2*m2 / mtot * np.dot(v1-v2, r1-r2) / d * (r1 - r2)
         u2 = v2 - 2*m1 / mtot * np.dot(v2-v1, r2-r1) / d * (r2 - r1)
         p[ii] = u1*m1
         p[jj] = u2*m2
           
     # Bounce off walls
     for i in range(Natoms+2):
        if pos[i,0] < Lmin:
            p[i,0] = -p[i,0]
            pos[i,0] = Lmin + abs(pos[i,0] % L)
        elif pos[i,0] > Lmax:
            p[i,0] = -p[i,0]
            pos[i,0] = Lmax - abs(pos[i,0] % L)
        if pos[i,1] < Lmin:
            p[i,1] = -p[i,1]
            pos[i,1] = Lmin + abs(pos[i,1] % L)
        elif pos[i,1] > Lmax:
            p[i,1] = -p[i,1]       
            pos[i,1] = Lmax - abs(pos[i,1] % L)
     
            
     tracex.append(pos[Natoms+1,0])
     tracey.append(pos[Natoms+1,1])
     
     t = t+dt
  
  # create fade line
  tracelc.remove()
  x = tracex[-4000:]
  y = tracey[-4000:]  
  lc = fade_line(x, y, colour, 1)
  tracelc = ax.add_collection(lc)      
  
  # statistics of travel distance
  x,y = circle(t)
  circleline.set_xdata(x)
  circleline.set_ydata(y)
   
  
  # Update plots  
  time = 0
  scatter1.set_xdata(pos[:Natoms,0])
  scatter1.set_ydata(pos[:Natoms,1]) 
  scatter2.set_xdata(pos[Natoms+1,0])
  scatter2.set_ydata(pos[Natoms+1,1]) 
  
   
  fig.suptitle('Time: ' + "{:.2f}".format(t/1e-3)+' ms')
    
  
  
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save('brownian.mp4',fps=25,dpi=300)
    
