#
# Lambert's law, emission from a surface
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
Natoms = 400  
L = 1 # size of volume in m
Matom = 4E-3/6E23 # helium mass
Ratom = 0.007 # size of helium atom in arbitrary units
k = 1.38E-23 
T = 300 # temperature
slitwidth = 0.08

dt = 5E-6
dtplot = 5e-6
simulationtime = 0.2
steps = int(simulationtime/dtplot)

NBbins = 50
thetascale = np.linspace(np.pi/2,-np.pi/2,NBbins)
thetabins = np.zeros(NBbins)

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
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(8,4.5))
fig.tight_layout(pad=0)

ax1 = plt.subplot(121)

slitx = [Lmax,Lmax]
slity = [(Lmax-Lmin)/2-slitwidth/2,(Lmax-Lmin)/2+slitwidth/2]

scatter1 = ax1.plot(pos[:,0],pos[:,1],marker='o',markersize=3,color='blue',lw=0)[0]
ax1.plot(slitx,slity,color='red',lw=5)
ax1.set(xlabel='position',ylabel='position',xlim=(Lmin,Lmax),ylim=(Lmin,Lmax))
ax1.set_aspect('equal')
text1 = ax1.title.set_text('Time:')

ax2 = plt.subplot(122, projection='polar')
ax2.set_thetamin(-90)
ax2.set_thetamax(90) 
ax2.set_yticks([])
ax2.set_ylim(0,1.1)
line1 = ax2.plot(thetascale,thetabins,color='blue')[0]
ax2.plot(thetascale,np.cos(thetascale),linestyle='dashed',color='r',label='cosine')
ax2.legend(fontsize=8)


t = 0.0
Nsteps = 0
time = 0
pos += (p/m)*(dt/2) # initial half-step
lostNB = 0

# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------

@jit(nopython = True)
def findhitlist(pos):
     hitlist = []
     # find collisions
     for i in range(Natoms+1):
        for j in range(i,Natoms+1):
           d2 = np.linalg.norm(pos[j]-pos[i])**2 
           if d2 < (2*Ratom)**2: 
               if i!=j:
                   hitlist.append((i,j)) 
     return hitlist              
     
def animate(k):
  global pos,p,t,time,m,lostNB

  print('step: '+ "{:.0f}".format(k)+' of '+"{:.0f}".format(steps)+', time:'+"{:.2f}".format(t/1e-3)) 
  
  #print('time: '+ "{:.2f}".format(t/1e-3)+' ms of '+"{:.0f}".format(simulationtime/1e-3)+' ms' )
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
           
     # Check for wall loss
     for i in range(Natoms+1):
        if pos[i,0] > Lmax and (pos[i,1] > (Lmax-Lmin)/2-slitwidth/2 and pos[i,1] < (Lmax-Lmin)/2+slitwidth/2):
           
           # sort partcile into distribution
           if p[i,0]>0:
             theta = np.arctan(p[i,1]/p[i,0])/(np.pi/2)*90
             print(theta)
             thetabini = int((90+theta)/180*NBbins)
             if thetabini >= 0 and thetabini < NBbins:
                 thetabins[thetabini] += 1
                 lostNB += 1
            
                 # create new particle
                 # pos[i,0] = Lmin+(Lmax-Lmin)*np.random.rand()
                 # pos[i,1] = Lmin+(Lmax-Lmin)*np.random.rand()
                
                 # r = Ratom
                 # mass = Matom*r**3/Ratom**3
                 # m = mass
                 # pavg = np.sqrt(2.*mass*1*k*T) # average kinetic energy p**2/(2mass) = (2/2)kT
    
                 # phi = 2*np.pi*np.random.rand()
                 # p[i,0] = pavg*np.cos(phi)
                 # p[i,1] = pavg*np.sin(phi)
        
     # Bounce off walls
     for i in range(Natoms+1):
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
  scatter1.set_xdata(pos[:,0])
  scatter1.set_ydata(pos[:,1])
  
  ma = np.max(thetabins)
  
  line1.set_ydata(thetabins/ma) 
    
  ax1.title.set_text('Time: ' + "{:.2f}".format(t/1e-3)+' ms, sampled particles: '+ "{:.0f}".format(lostNB))
    
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save('lambert.mp4',fps=25,dpi=300)
    
