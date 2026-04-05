#
# Diffusive mixing of two gases
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#

from numba import jit
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

# ------------------------------------------------------------------------
# parameter
# ------------------------------------------------------------------------ 

Natoms = 300  
L = 1 # size of volume in m
Matom = 4E-3/6E23 # helium mass
Ratom = 0.007 # size of helium atom in arbitrary units
k = 1.38E-23 
T = 300 # temperature

dt = 5E-6
dtplot = 5e-6
simulationtime = 1e-2
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
    x = Lmin+(Lmax-Lmin)/2*np.random.rand()
    y = Lmin+(Lmax-Lmin)*np.random.rand()
    r = Ratom
    mass = Matom*r**3/Ratom**3
    m = mass
    pavg = np.sqrt(2.*mass*1*k*T) # average kinetic energy p**2/(2mass) = (2/2)kT
    
    phi = 2*np.pi*np.random.rand()
    px = pavg*np.cos(phi)
    py = pavg*np.sin(phi)

    poslist.append((x,y))
    plist.append((px,py))
    mlist.append(mass)
    rlist.append(r)
    
for i in range(Natoms+1):
    x = L/2+(Lmax-Lmin)/2*np.random.rand()
    y = Lmin+(Lmax-Lmin)*np.random.rand()
    r = Ratom
    mass = Matom*r**3/Ratom**3
    m = mass
    pavg = np.sqrt(2.*mass*1*k*T) # average kinetic energy p**2/(2mass) = (2/2)kT
    
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
fig, ax = plt.subplots(1,2,figsize=(8,4.5))
fig.tight_layout(pad=3)

scatter1 = ax[0].plot(pos[:Natoms,0],pos[:Natoms,1],marker='o',markersize=3,color='blue',lw=0)[0]
scatter2 = ax[0].plot(pos[Natoms:,0],pos[Natoms:,1],marker='o',markersize=3,color='red',lw=0)[0]
ax[0].set(xlabel='position',ylabel='position',xlim=(Lmin,Lmax),ylim=(Lmin,Lmax))
ax[0].set_aspect('equal')


line1 = ax[1].plot(xscale,xbin1,color='blue')[0]
line2 = ax[1].plot(xscale,xbin2,color='red')[0]
ax[1].set(xlabel='position',ylabel='N(x) (norm)',ylim=(0,0.1),xlim=(0,L))
ax[1].legend()


t = 0.0
Nsteps = 0
time = 0
pos += (p/m)*(dt/2) # initial half-step


# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------

@jit(nopython = True)
def findhitlist(pos):
     hitlist = []
     # find collisions
     for i in range(2*Natoms+1):
        for j in range(i,2*Natoms+1):
           d2 = np.linalg.norm(pos[j]-pos[i])**2 
           if d2 < (2*Ratom)**2: 
               if i!=j:
                   hitlist.append((i,j)) 
     return hitlist              
     
def animate(k):
  global pos,p,t,time

  print('time: '+ "{:.2f}".format(t/1e-6)+' microseconds of '+"{:.0f}".format(simulationtime/1e-6)+' microseconds' )
  t += dt
  while time<=dtplot:
     time += dt   

     pos += p/m*dt   
     
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
         u1 = v1 - 2*m1 / mtot * np.dot(v1-v2, r1-r2) / d * (r1 - r2)
         u2 = v2 - 2*m2 / mtot * np.dot(v2-v1, r2-r1) / d * (r2 - r1)
         p[ii] = u1*m1
         p[jj] = u2*m2
           
     # Bounce off walls
     for i in range(2*Natoms+1):
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
               
     t = t+dt
  
  # Update plots  
  time = 0
  scatter1.set_xdata(pos[:Natoms,0])
  scatter1.set_ydata(pos[:Natoms,1]) 
  scatter2.set_xdata(pos[Natoms:,0])
  scatter2.set_ydata(pos[Natoms:,1]) 
  
  for i in range(NBbins):
      xbin1[i] = 0        
  for i in range(Natoms):
      xi = pos[i,0]
      binx = round(xi/(L/NBbins)) 
      if binx < NBbins:
          xbin1[binx] += 1/Natoms 
  line1.set_ydata(xbin1) 
  
  for i in range(NBbins):
      xbin2[i] = 0        
  for i in range(Natoms):
      xi = pos[Natoms+i,0]
      binx = round(xi/(L/NBbins)) 
      if binx < NBbins:
          xbin2[binx] += 1/Natoms 
  line2.set_ydata(xbin2) 
  
    
  fig.suptitle('Time: ' + "{:.2f}".format(t/1e-3)+' ms')
    
  
  
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save('diffusion.mp4',fps=25,dpi=300)
    
