#
# Thermalization of a distribution function
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#


from numba import jit
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

# ------------------------------------------------------------------------
# parameter
# ------------------------------------------------------------------------ 
scaling = 'plotvsE'  # plot distribution vs. E
scaling = 'plotvsv' # plot distribution vs. |v|

if scaling == 'plotvsv':
  animationname = 'gaskineticvsv.mp4'
if scaling == 'plotvsE':
  animationname = 'gaskineticvsE.mp4'

Natoms = 1000  
L = 1 # size of volume in m
Matom = 4E-3/6E23 # helium mass
Ratom = 0.01 # size of helium atom in arbitrary units
k = 1.38E-23 
T = 300 # temperature

dt = 5E-6
dtplot = 5e-6
simulationtime = 2e-3
steps = int(simulationtime/dtplot)

NBbins = int(Natoms/5)
if scaling == 'plotvsE':
    Ebin = np.zeros(NBbins)
    ERange = 0.25*1e-19
    binwidth = ERange/NBbins
    Escale = np.linspace(0,ERange,NBbins) 
elif scaling == 'plotvsv':    
    vbin = np.zeros(NBbins)
    vRange = np.sqrt(0.25*1e-19*2/Matom)
    binwidth = vRange/NBbins
    vscale = np.linspace(0,vRange,NBbins) 


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
    
    phi = 0.1*np.pi*np.random.rand()
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

scatter1 = ax[0].plot(pos[:,0],pos[:,1],marker='o',markersize=3,lw=0)[0]
ax[0].set(xlabel='position',ylabel='position',xlim=(Lmin,Lmax),ylim=(Lmin,Lmax))
ax[0].set_aspect('equal')


if scaling == 'plotvsE':
    nemax = 0.06
    line1 = ax[1].plot(Escale/1e-19,Ebin)[0]
    ax[1].set(xlabel='energy (eV)',ylabel='N(E) (norm)',ylim=(0,nemax),xlim=(0,NBbins*binwidth/1e-19))

    Etheory = np.linspace(0,ERange,1000)
    NEtheory = np.zeros(1000)
    for i in range(1000):
        vt = np.sqrt(2*Etheory[i]/(m))
        NEtheory[i] = 1/(k*T)*np.exp((-0.5*m*vt**2)/(k*T))*binwidth
    line2 = ax[1].plot(Etheory/1e-19,NEtheory,color='r',linestyle='dashed',label='theory')[0]
    ax[1].legend()
elif scaling == 'plotvsv':
    nemax = 0.03    
    line1 = ax[1].plot(vscale,vbin)[0]
    ax[1].set(xlabel='velocity (m/s)',ylabel='N(v) (norm)',ylim=(0,nemax),xlim=(0,NBbins*binwidth))

    vtheory = np.linspace(0,vRange,1000)
    nvtheory = np.zeros(1000)
    for i in range(1000):
        vt = vtheory[i]
        nvtheory[i] = m/(k*T)*np.exp((-0.5*m*vt**2)/(k*T))*vt*binwidth
    line2 = ax[1].plot(vtheory,nvtheory,color='r',linestyle='dashed',label='theory')[0]
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
     for i in range(Natoms+1):
        for j in range(i,Natoms+1):
           d2 = np.linalg.norm(pos[j]-pos[i])**2 
           if d2 < (2*Ratom)**2: 
               if i!=j:
                   hitlist.append((i,j)) 
     return hitlist     

def animate(k):
  global pos,p,t,time

  print('time: '+ "{:.2f}".format(t/1e-6)+' microseconds of '+"{:.0f}".format(simulationtime/1e-6)+' microseconds' )
  t += dt
  while time<=dtplot-dt:
     time += dt   

     pos += p/m*dt   
     
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

     hitlist =  findhitlist(pos)
     #print(hitlist)
     hits = np.array(hitlist)           
  
     # change momentum
     for i in range(len(hits)):
         ii = hits[i,0]
         jj = hits[i,1]         
      #   print(i)    
         mtot = m+m    
    
         r1 = pos[ii]
         r2 = pos[jj]
         v1 = p[ii]/m
         v2 = p[jj]/m
         d = np.linalg.norm(r1 - r2)**2
         u1 = v1 - 2*m / mtot * np.dot(v1-v2, r1-r2) / d * (r1 - r2)
         u2 = v2 - 2*m / mtot * np.dot(v2-v1, r2-r1) / d * (r2 - r1)
         p[ii] = u1*m
         p[jj] = u2*m
           

               
     t = t+dt
  
  # Update plots  
  time = 0
  scatter1.set_xdata(pos[:,0])
  scatter1.set_ydata(pos[:,1]) 
  
  if scaling == 'plotvsE':
      for i in range(NBbins):
          Ebin[i] = 0
      for i in range(Natoms):
          ei = np.linalg.norm(p[i])**2/(2*m)
          binenergy = round(ei/binwidth) 
          if binenergy < NBbins:
              Ebin[int(ei/binwidth)] += 1 
      for i in range(NBbins):
       if Ebin[i]/Natoms>nemax:
            Ebin[i]=nemax*Natoms 
      line1.set_xdata(Escale/1e-19)
      line1.set_ydata(Ebin/Natoms) 
  if scaling == 'plotvsv':
      for i in range(NBbins):
          vbin[i] = 0
      for i in range(Natoms):
          vi = np.linalg.norm(p[i])/m
          binv = round(vi/binwidth) 
          if binv < NBbins:
              vbin[int(vi/binwidth)] += 1 
      for i in range(NBbins):
       if vbin[i]/Natoms>nemax:
            vbin[i]=nemax*Natoms 
      line1.set_xdata(vscale)
      line1.set_ydata(vbin/Natoms) 
      
  fig.suptitle('Time: ' + "{:.0f}".format(t/1e-6)+' microseconds')
    
  
  
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)
    
