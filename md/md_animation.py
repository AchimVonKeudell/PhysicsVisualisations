#
# Molceular Dynamic Simulation Metal
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

animationname = 'md.mp4'

el = 1.6e-19
amu = 1.6e-27
kB = 1.38e-23
t = 0
time = 0

tend = 1e-12
deltat = 5e-16
dtreplot = 1e-15
steps = int(tend/dtreplot)
ProjectileEnergy = 500 # energy in eV
masstarget = 26
massprojectile = 26
temperature = 300 # temperature bath outer atoms
tautemp = 2e-15 # energy relaxation time

BasisNB = 2
xNB = 12
yNB = 12
zNB = 8
maxnb = xNB*yNB*zNB*BasisNB+1

dneighbours = 5*1e-10
RPotCut = 10*1e-10
rskin = 2*1e-10

class atom:
    def __init__(self,r,dr,v,a,Epot,Ekin,name,mass,neighbour,neighbourNB):
        self.r = r # location
        self.dr = dr # 
        self.v = v # velocity
        self.a = a # acceleration
        self.Epot = Epot
        self.Ekin = Ekin
        self.name = name
        self.mass = mass
        self.neighbour = neighbour
        self.neighbourNB = neighbourNB
  
atoms = []
for i in range(maxnb):
  atoms.append(atom(r=np.zeros(3),dr=np.zeros(3),v=np.zeros(3),a=np.zeros(3),
                   Ekin=0,Epot=0,name='atom',mass=1,
                   neighbour=np.zeros(50),neighbourNB=10))

class basisatom:
    def __init__(self,r):
        self.r = r # location
          
basis = []
for i in range(8):
  basis.append(basisatom(r=np.zeros(3)))


  
def VectorPeriodic(v):
  global xsize,ysize,zsize  
  if v[0] < -(xsize/2): v[0] += + xsize
  if v[0] >= (xsize/2): v[0] -= - xsize
  if v[1] < -(ysize/2): v[1] += + ysize
  if v[1] >= (ysize/2): v[1] -= - ysize
  #if v[2] < -(zsize/2): v[2] += + zsize
  #if v[2] >= (zsize/2): v[2] -= - zsize
  return v

# ------------------------------------------------------------------------
# Find neighbours
# ------------------------------------------------------------------------
def GenerateNeighbours():
  global AtomNB  
  for i in range(AtomNB):
      atoms[i].NeighbourNB = 0
      for j in range(AtomNB):
         if i!=j:
              dr = atoms[j].r-atoms[i].r
              dr = VectorPeriodic(dr)
              drb = np.linalg.norm(dr);
              if (drb!=0) and (drb < dneighbours):
                  atoms[i].neighbour[atoms[i].NeighbourNB] = j
                  atoms[i].NeighbourNB += 1


# --------------------------------------------------------------------------
# initialize atoms
# -------------------------------------------------------------------------

# target
alattice = 3.75*1e-10   
xsize = (xNB)*alattice
ysize = (yNB)*alattice
zsize = (zNB)*alattice

# displayvolume in AA
xylim = xsize/1e-10*0.5
zlim = 2*zsize/1e-10*0.5
 
# Projectile
ProjectileV = np.sqrt(ProjectileEnergy*el*2/(28*amu));

atoms[0].r = np.array([1e-10,2e-10,zlim*1e-10])
atoms[0].dr = np.zeros(3)
atoms[0].v = np.array([0,0,-ProjectileV]) 
atoms[0].a = np.zeros(3)
atoms[0].name = 'pro'
atoms[0].mass = massprojectile*amu
    
# generate basis bcc
basis[0].r = np.zeros(3)
basis[1].r = np.array([0.5,0.5,0.5])
BasisNB = 2
  
# target
vrms = np.sqrt(kB*temperature/(masstarget*amu))
vrms = 0
  
i = 0
corner = [-xsize/2,-ysize/2,-zsize/2]
for x in range(xNB):
    for y in range(yNB):
         for z in range(zNB):
            for k in range(BasisNB):
                   i += 1
                   offset = np.array([x,y,z])
                   offset2 = np.array([0,0,0])
                   atoms[i].r = corner + (offset+basis[k].r+offset2)*alattice                   
                   atoms[i].dr = np.array([0,0,0])
                   atoms[i].v = vrms*np.array([np.random.normal(),np.random.normal(),np.random.normal()])                 
                   atoms[i].a = np.array([0,0,0])
                   atoms[i].name = 'Al'
                   atoms[i].mass = masstarget*amu

AtomNB = i+1
GenerateNeighbours()
                

# -------------------------------------------------------------------
# Morse Potential
# -------------------------------------------------------------------
MorseDe = 0.2703*el
MorseAlpha = 1.1646*1/1e-10
MorseRe = 3.253*1e-10

def Morse(x):
   ex = np.exp(-MorseAlpha*(x-MorseRe))  
   return -MorseDe*(1-ex)**2
  
def Morsedr(x):
  ex = np.exp(-MorseAlpha*(x-MorseRe))
  return -2*MorseDe*MorseAlpha*(ex**2-ex)


# ---------------------------------------------------------------------------
# Generate Forces
# ---------------------------------------------------------------------------
def GenerateForces2body():
  global AtomNB  
  #print(AtomNB)
  for i in range(AtomNB):
      atoms[i].a = np.array([0,0,0])
      atoms[i].Epot = 0
  
  for i in range(AtomNB):
      AtomPt1 = i
      for j in range(atoms[i].NeighbourNB):
          AtomPt2 = int(atoms[i].neighbour[j])
          rij = atoms[AtomPt2].r-atoms[AtomPt1].r
          rij = VectorPeriodic(rij)
          rijb = np.linalg.norm(rij)
          #print(rijb)
          if (rijb<RPotCut):
              V    = Morse(rijb)
              dVdr = Morsedr(rijb)
              a    = -dVdr/atoms[AtomPt1].mass
              #print(a)
              atoms[AtomPt1].a = atoms[AtomPt1].a - rij/rijb*a;
              atoms[AtomPt1].Epot += V
              atoms[AtomPt2].a = atoms[AtomPt2].a + rij/rijb*a;


# ------------------------------------------------------------------------
# Setup plot
# ------------------------------------------------------------------------
xt = []
yt = []
zt = []
fig = plt.figure(figsize=(8,4.5))
axp = fig.add_subplot(1,1,1,projection='3d')
axp.set(xlabel="x (AA)",ylabel="y (AA)",zlabel="z (AA)")
axp.set(xlim=(-xylim,xylim),ylim=(-xylim,xylim),zlim=(-zlim,zlim))
axp.set_box_aspect((1,1,zlim/xylim))
#axp.quiver(0,0,0,0,0,1,pivot = 'tip', color="r", length=20,label='impact',alpha=0.5)

# trajectory line
linet = axp.plot3D(xt,yt,zt,color='b',marker='o',markersize=3,label='target',lw=0)[0]
linep = axp.plot3D(xt,yt,zt,color='r',marker='o',markersize=3,label='projectile',lw=0)[0]
axp.legend(fontsize=6)


# info box
infobox = ''
infobox += 'E proj.: ' + "{:.0f}".format(ProjectileEnergy) + ' (eV)\n' 
infobox += 'm proj.: ' + "{:.0f}".format(massprojectile) + ' (amu)\n' 
infobox += 'm target: ' + "{:.0f}".format(masstarget) + ' (amu)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
axp.text2D(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=axp.transAxes)


# ---------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------

def animate(k):

  global time, t
  
  # Search for neighbours  
  GenerateNeighbours()
  
  # propagation step 1
  time = 0
  while time<dtreplot:
    t += deltat
    time += deltat  
    for i in range(AtomNB):
      dr = atoms[i].v*deltat+0.5*deltat**2*atoms[i].a
      #if i>0: print(dr)

      atoms[i].dr = atoms[i].dr + dr
      atoms[i].r = atoms[i].r + dr
      atoms[i].v = atoms[i].v + 0.5 * deltat*atoms[i].a

      #atoms[i].r = VectorPeriodic(atoms[i].r)

    # Generate forces
    GenerateForces2body()
  
    # second half kick
    for i in range(AtomNB):
      atoms[i].v = atoms[i].v + 0.5*deltat*atoms[i].a
      atoms[i].Ekin = 0.5*3.75e2*atoms[i].mass*(atoms[i].v[0]**2+atoms[i].v[1]**2+atoms[i].v[2]**2)

    # find energy outer atoms 
    ekinsum = 0
    skinnb = 0
    for i in range(AtomNB):
        if (atoms[i].r[0]<(-xsize/2+rskin) or atoms[i].r[0]>(xsize/2-rskin) or
            atoms[i].r[1]<(-ysize/2+rskin) or atoms[i].r[1]>(ysize/2-rskin)):
           skinnb += 1
           ekinsum += atoms[i].v[0]**2+atoms[i].v[1]**2+atoms[i].v[2]**2    
    tempskin = 0.5*atoms[1].mass*ekinsum/(1.5*AtomNB*kB)
    tempscale = np.sqrt(1+deltat/tautemp*(temperature/tempskin-1));
    
    # rescale outer atoms 
    for i in range(AtomNB):
        if (atoms[i].r[0]<(-xsize/2+rskin) or atoms[i].r[0]>(xsize/2-rskin) or
            atoms[i].r[1]<(-ysize/2+rskin) or atoms[i].r[1]>(ysize/2-rskin)):
           atoms[i].v = atoms[i].v * tempscale    
    
  # ------------------------------------------------------
  # Display solution
  # ------------------------------------------------------    
  # target
  x1 = []
  y1 = []
  z1 = []
  for i in range(1,AtomNB):
        x1.append(atoms[i].r[0]/1e-10)
        y1.append(atoms[i].r[1]/1e-10)
        z1.append(atoms[i].r[2]/1e-10)
            
  linet.set_xdata(x1)
  linet.set_ydata(y1) 
  linet.set_3d_properties(z1)
  
  # projectile
  linep.set_xdata([atoms[0].r[0]/1e-10])
  linep.set_ydata([atoms[0].r[1]/1e-10]) 
  linep.set_3d_properties([atoms[0].r[2]/1e-10])
          
  axp.view_init(20,35+t/tend*180)        
  
  # Update the title of the plot
  fig.suptitle('time: ' + "{:.0f}".format(t/1e-15)+' fs')
  print('time: ' + "{:.0f}".format(t/1e-15)+' fs of '+ "{:.0f}".format(tend/1e-15)+' fs')

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=180)
 
