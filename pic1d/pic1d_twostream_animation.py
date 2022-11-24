#
# Two Stream instabilitiy in a 1d1v PIC code
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#model = 'normalizedunits'
model = 'SIunits'

NBParticles = 70000
NBBins = 1000

if model == 'normalizedunits':
   xmin = 0
   xmax = 100 
   dt = 0.05
   dtplot = 0.1
   tend = 3
   Escale = 5
   vscale = 6
   el = 1
   me = 1
   eps0 = 1
   vth = 1
   qsuper = 1
   n0 = 1
   kB = 1.38e-23
if model == 'SIunits':
   xmin = 0
   xmax = 0.01
   dt = 1e-10
   dtplot = 2e-10
   tend = 5e-9

   el = 1.6e-19
   me = 6e-31
   kB = 1.38e-23
   n0 = 1e15
   eps0 = 8e-12
   Te = 3
   vth =  0.02*np.sqrt(3*el*Te/(me))
   vscale = 6*vth/1e3    
   Escale = 100
   qsuper = n0*NBBins/NBParticles

ldebye = np.sqrt(eps0*kB*Te/(n0*el**2))
fpe = np.sqrt(n0*el**2/(eps0*me))*2*np.pi

deltax = (xmax-xmin)/NBBins
widthx = NBBins*deltax
averagen0 = 1

class picbin:
    def __init__(self,n,x,E):
        self.n = n
        self.x = x
        self.E = E
  
bins = []
for i in range(NBBins):
  bins.append(picbin(n=0,x=0,E=0))

class picparticle:
    def __init__(self,x,v,a,E,color):
        self.x = x
        self.v = v
        self.a = a
        self.E = E
        self.color = color

particles = []
for i in range(NBParticles):
  particles.append(picparticle(x=0,v=0,a=0,E=0,color=0))                                               

def randomnormal():
 return np.sqrt(-2*np.log(np.random.rand()))*np.cos(2*3.14159*np.random.rand());

def ParticlesToBins():
  for j in range(NBBins):
      bins[j].n = 0
  totaln = 0

  for i in range(NBParticles):
     jbin = int(np.trunc(particles[i].x/deltax))-1
     jbinp1 = jbin + 1
     #print(jbinp1)
     bins[jbin].n += 1/deltax*(bins[jbinp1].x-particles[i].x)
     bins[jbinp1].n += 1/deltax*(particles[i].x-bins[jbin].x)
     totaln += 1/deltax*(bins[jbinp1].x-particles[i].x)+1/deltax*(particles[i].x-bins[jbin].x)
   
  for j in range(NBBins):
      bins[j].n = bins[j].n*NBBins/NBParticles*averagen0-averagen0

def CalculateEfield():
  bins[NBBins-1].E =0 
  Esum = 0 
  for j in range(NBBins-1,1,-1):
      Esum += bins[j].E;
      bins[j-1].E = bins[j].E + 0.5*deltax*bins[j].n*qsuper*el/eps0+0.5*deltax*bins[j-1].n*qsuper*el/eps0

def EFieldToParticles():
  for i in range(NBParticles):
      jbin = int(np.trunc(particles[i].x/deltax))-1
      jbinp1 = jbin+1;
      particles[i].E = bins[jbin].E*(bins[jbinp1].x-particles[i].x)/deltax+bins[jbinp1].E*(particles[i].x-bins[jbin].x)/deltax

def ParticlesAcceleration():
  for i in range(NBParticles):
      particles[i].a = -particles[i].E*el/me;

def InitializeParticles():
  for i in range(NBBins):
      bins[i].x = 0.5*deltax+(i-1)*deltax
      bins[i].n = 0
  vmax = 1

  for i in range(NBParticles):
      particles[i].x = np.random.rand()*widthx;
      particles[i].v = vth*(2+1*randomnormal())*(1 + 0.1*np.sin(6*3.14159*(particles[i].x)/widthx))
      particles[i].a = 0
      particles[i].color = 1
   
  for i in range(int(NBParticles*0.5),NBParticles):
      particles[i].v = particles[i].v*(-1)
      particles[i].color = 2

InitializeParticles()
ParticlesToBins()

# --------------------------------------------------------------------------
# Define plots
# --------------------------------------------------------------------------
fig, ax = plt.subplots(1,1,figsize=(8,4.5))
plt.subplots_adjust(left=0.2, right=0.8, top=0.85, bottom=0.1)


x1 = []
for i in range(int(NBParticles/2)):
    x1.append(particles[i].x)
v1 = []
for i in range(int(NBParticles/2)):
    v1.append(particles[i].v)
x2 = []
for i in range(int(NBParticles/2),NBParticles):
    x2.append(particles[i].x)
v2 = []
for i in range(int(NBParticles/2),NBParticles):
    v2.append(particles[i].v)
xe = []
for i in range(NBBins):
    xe.append(bins[i].x)
E = []
for i in range(NBBins):
    E.append(bins[i].E)


ax.set(xlabel='x (m)',ylabel='v (km/s)',ylim=(-vscale,vscale),xlim=(xmin,xmax))
line1 = ax.plot(x1,v1,label='ne',color='r',markersize=0.5,marker='o',lw=0)[0]
line2 = ax.plot(x2,v2,label='ne',color='b',markersize=0.5,marker='o',lw=0)[0]
axp = ax.twinx()
axp.set(ylabel='E field (V/m)',ylim=(-Escale,Escale))
lineE = axp.plot(xe,E,label='E',color='g',lw=2,alpha=0.8)[0]
ax.legend(loc=3,fontsize=6)
axp.legend(loc=1,fontsize=6)

# info box
infobox = ''
infobox += 'n0: ' + "{:.0e}".format(n0) + ' (1/m^3)\n'
infobox += 'l_Debye: ' + "{:.0e}".format(ldebye) + ' (m)\n'    
infobox += 'fpe: ' + "{:.0e}".format(fpe) + ' (1/s)'    
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8) 
ax.text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=axp.transAxes)



t = 0
steps = int(tend/dtplot)
animationname = '2streaminstability.mp4'

# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
  global t,v,x,a
  global xt,yt,zt
  
  time = 0
  while time < dtplot:    
    t += dt
    time += dt
    
    # 1st half kick 
    for i in range(NBParticles):
        particles[i].v += particles[i].a*0.5*dt

    # Update positions 
    for i in range(NBParticles):
        particles[i].x += particles[i].v*dt;
        # Periodic boundary condisitons 
        if particles[i].x>xmax:
            particles[i].x = particles[i].x-xmax;
        if particles[i].x<xmin:
            particles[i].x = xmax-abs(particles[i].x)

    # new Accelerations 
    ParticlesToBins()
    CalculateEfield()
    EFieldToParticles()
    ParticlesAcceleration()

    # second half kick 
    for i in range(NBParticles):
        particles[i].v += particles[i].a*0.5*dt;
     
    
  x1 = []
  for i in range(int(NBParticles/2)):
    x1.append(particles[i].x)
  v1 = []
  for i in range(int(NBParticles/2)):
    v1.append(particles[i].v/1e3)
  x2 = []
  for i in range(int(NBParticles/2),NBParticles):
    x2.append(particles[i].x)
  v2 = []
  for i in range(int(NBParticles/2),NBParticles):
    v2.append(particles[i].v/1e3)
  E = []
  for i in range(NBBins):
    E.append(bins[i].E)

    
  line1.set_xdata(x1)
  line1.set_ydata(v1) 
  line2.set_xdata(x2)
  line2.set_ydata(v2)
  lineE.set_ydata(E)
  #z = []
  #for i in range(NBBins):
  #  z.append([bins[i].E,0])
  #cp.set_3dproperties(z)
        
  # Update the title of the plot
  fig.suptitle('time: ' + "{:.1f}".format(t/1e-9)+' ns')
  print('time: ' + "{:.1f}".format(t/1e-9)+' ns of '+ "{:.1f}".format(tend/1e-9)+' ns')
           
       
# ----------------------------------------
# Create the complete time sequence
# the total time span is tend/dtplot
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=180)
                