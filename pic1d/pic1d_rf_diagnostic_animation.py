#
# plasma in a 1d1v PIC code
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec

#model = 'plasmaonset'
#model = 'dc'
model = 'rf'
#model = 'doublelayer'

NBParticles = 10000
NBBins = 200

NBEEDF = 100
EEDFmax = 15


if model == 'plasmaonset':
    Te = 3
    V0 = 0
    E0y = 0
    phiscale = 10*Te
    dt = 1e-10
    dtplot = 2e-10
    tend = 500e-9
    n0 = 1e15
    nuel = 1e7
    nuion = 1e6
    animationname = 'plasmaonset.gif'
    useinitial = True
    initialfilename = 'plasmaonset.dat'
if model == 'dc':
    Te = 3
    V0 = -100
    E0y = 0
    phiscale = 12*Te
    dt = 1e-10
    dtplot = 2e-10
    tend = 500e-9
    n0 = 1e15
    nuel = 1e7
    nuion = 1e6
    animationname = 'dcm100Vnu1e7el1e6ion.mp4'
    useinitial = False
    initialfilename = 'dcm100.dat'
if model == 'rf':
    Te = 3
    V0 = -100
    E0y = 0
    phiscale = 12*Te
    dt = 1e-10
    dtplot = 2e-10
    tend = 5*74e-9 # 2e-10#
    n0 = 1e15
    nuel = 1e7
    nuion = 0#1e6
    animationname = 'rfplasma_diag.mp4'
    useinitial = False
    initialfilename = 'rfplasma.dat'
if model == 'doublelayer':
    Te = 3
    V0 = 0
    E0y = 20000
    phiscale = 30*Te
    dt = 1e-10
    dtplot = 2e-10
    tend = 500e-9
    n0 = 1e15
    nuel = 1e7
    nuion = 1e6
    animationname = 'doublelayer.mp4'
    useinitial = False
    initialfilename = 'doublelayer.dat'

xmin = 0
xmax = 0.025

el = 1.6e-19
me = 6e-31
mion = 1.6e-27*40
massratio1 = me/(me+mion)
massratio2 = mion/(me+mion)
kB = 1.38e-23
eps0 = 8.854e-12
Tg = 300
frf = 13.45e6
vthelectrons =  np.sqrt(3*el*Te/(me))
vthions =  0#np.sqrt(3*kB*Tg/(mion))
vscale = 15*vthelectrons    


qsuper = n0/NBParticles*NBBins

print('Start simulation')

ldebye = np.sqrt(eps0*kB*Te/(n0*el**2))
fpe = np.sqrt(n0*el**2/(eps0*me))*2*np.pi

deltax = (xmax-xmin)/(NBBins)
widthx = NBBins*deltax

class picbin:
    def __init__(self,rho,npos,nneg,x,E,phi):
        self.rho = rho # net charge
        self.npos = npos # positive charge
        self.nneg = nneg # negative charge
        self.x = x
        self.E = E
        self.phi = phi
  
bins = []
for i in range(NBBins):
  bins.append(picbin(rho=0,npos=0,nneg=0,x=0,E=0,phi=0))
f = np.zeros(NBBins)
w = np.zeros(NBBins)
g = np.zeros(NBBins)

class picparticle:
    def __init__(self,x,vx,vy,vz,ax,ay,E,q,m,label):
        self.x = x
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.ax = ax
        self.ay= ay
        self.E = E
        self.q = q
        self.m = m
        self.label = label

particles = []
for i in range(NBParticles):
  particles.append(picparticle(x=0,vx=0,vy=0,vz=0,ax=0,ay=0,E=0,q=1,m=1,label=0))                                               

def voltageelectrode(t):
    if model == 'plasmaonset':
        return V0
    if model == 'dc':
        return V0
    if model == 'rf':
        return V0*np.cos(2*np.pi*frf*t) 
    if model == 'doublelayer':
        return V0

def Ey(t):
    return E0y*np.cos(2*np.pi*frf*t) 
    
def randomnormal():
 return np.sqrt(-2*np.log(np.random.rand()))*np.cos(2*3.14159*np.random.rand());

def ParticlesToBins():
  for j in range(NBBins):
      bins[j].rho = 0
      bins[j].npos = 0
      bins[j].nneg = 0

  for i in range(NBParticles):
     jbin = int(np.trunc(particles[i].x/deltax))
     if jbin<NBBins-1:
       jbinp1 = jbin + 1
       bins[jbin].rho += 1/deltax*(bins[jbinp1].x-particles[i].x)*particles[i].q
       if particles[i].label == 2:
           bins[jbin].npos += 1
       elif particles[i].label == 1:
           bins[jbin].nneg += 1 
       bins[jbinp1].rho += 1/deltax*(particles[i].x-bins[jbin].x)*particles[i].q
       if particles[i].label == 2:
           bins[jbinp1].npos += 1
       elif particles[i].label == 1:
           bins[jbinp1].nneg += 1
     else:
       bins[jbin].rho += 1/deltax*(xmax-particles[i].x)*particles[i].q
         

# --------------------------------------------------------------------------
# Solve Poissons equation by direct inversion of matrix
# code from eduPIC
# --------------------------------------------------------------------------     
def CalculateEfield2(t):
  global bins  
  A = 1.0
  B = -2.0;
  C = 1.0
  S = 1.0 / (2.0 * deltax);
  ALPHA = -deltax**2 / eps0

  bins[0].phi = voltageelectrode(t)#V0 #VOLTAGE * cos ( OMEGA * tt); // potential at the powered electrode
  bins[NBBins -1].phi = 0 #; // potential at the grounded electrode

  # solve Poisson equation
 
  for i in range(1,NBBins):
    f[i] = ALPHA * bins[i].rho*qsuper*el;
  f[1] -= bins[0].phi;
  f[NBBins -2] -= bins[NBBins -1].phi
  w[1] = C/B;
  g[1] = f [1]/ B;
  for i in range(2,NBBins-1):
    w[i] = C / (B - A * w[i -1]) ;
    g[i] = (f[i] - A * g[i -1]) / (B - A * w[i -1]) ;

  bins[NBBins -2].phi = g[NBBins -2];
  for i in range(NBBins -2,0,-1):
      bins[i].phi = g[i] - w[i] * bins[i +1].phi #// potential at the grid points between the electrodes

  for i in range(1,NBBins-1):
      bins[i].E = (bins[i -1].phi - bins[i +1].phi) * S #; // electric field at the grid points between the electrodes

  bins[0].E = (bins[0].phi - bins[1].phi) * 1/deltax -  bins[0].rho*qsuper*el * deltax / (2.0 * eps0 ); # powered electrode
  bins[NBBins -1].E = ( bins[NBBins -2].phi - bins[NBBins -1].phi) * 1/deltax  +  bins[NBBins-1].rho*qsuper*el * deltax / (2.0 * eps0); # grounded
  

def EFieldToParticles():
  for i in range(NBParticles):
      jbin = int(np.trunc(particles[i].x/deltax))
      if jbin<NBBins-1:
        jbinp1 = jbin+1;
        particles[i].E = bins[jbin].E*(bins[jbinp1].x-particles[i].x)/deltax+bins[jbinp1].E*(particles[i].x-bins[jbin].x)/deltax
      else:  
        particles[i].E = bins[jbin].E#*(xmax-particles[i].x)/deltax#+bins[jbinp1].E*(particles[i].x-bins[jbin].x)/deltax
    
def ParticlesAcceleration(t):
  # acceleration in 1d direction  
  for i in range(NBParticles):
      particles[i].ax = particles[i].E*particles[i].q*el/particles[i].m;
  # loacl rf heating in y direction   
  for i in range(NBParticles):    
      if particles[i].label == 1 and particles[i].x<xmax*0.5 and particles[i].x>xmax*0.4: 
          particles[i].ay = Ey(t)*particles[i].q*el/particles[i].m;


def InitializeParticles():
  global NBParticles
  
  for i in range(NBBins):
      bins[i].x = 0.5*deltax+i*deltax
      bins[i].n = 0
  
  # electrons  
  for i in range(NBParticles):
      particles[i].x = np.random.rand()*widthx;
      particles[i].vx = vthelectrons*randomnormal()
      #particles[i].vx = np.random.normal(scale=np.sqrt(3*el*Te/me))      
      particles[i].vy = 0
      particles[i].vz = 0
      particles[i].ax = 0
      particles[i].q = -1
      particles[i].m = me
      particles[i].label = 1
   
  # ions   
  for i in range(int(NBParticles*0.5),NBParticles):
      particles[i].vx = vthions*randomnormal()
      particles[i].q = 1
      particles[i].m = mion
      particles[i].label = 2
  
  if useinitial:   
      data = np.loadtxt(initialfilename)     
      # ,[xp,vxp,vyp,qp,mp,labelp]
      NBParticles = len(data[1])
      for i in range(NBParticles):
          particles[i].x = data[0,i]
          particles[i].vx = data[1,i]
          particles[i].vy = data[2,i]
          particles[i].q = data[3,i]
          particles[i].m = data[4,i]
          particles[i].label = data[5,i]

# --------------------------------------------------------------------------
# Begin setting up simulation
# --------------------------------------------------------------------------
InitializeParticles()
ParticlesToBins()

# --------------------------------------------------------------------------
# Define plots
# --------------------------------------------------------------------------
#fig, ax = plt.subplots(2,1,figsize=(8,4.5),gridspec_kw={'height_ratios': [2,1]})
fig = plt.figure(constrained_layout=True)
ax = []
gs = GridSpec(2,2,figure=fig,height_ratios = (2,1),width_ratios = (1,1))
ax.append(fig.add_subplot(gs[0,0])) # phase space
ax.append(fig.add_subplot(gs[1,0])) # density profile
ax.append(fig.add_subplot(gs[0:,-1])) # EEDF
ax[2].yaxis.tick_right()

x1 = []
v1 = []
x2 = []
v2 = []
for i in range(NBParticles):
    if particles[i].label == 1: # electrons
        x1.append(particles[i].x/1e-3)
        v1.append(particles[i].vx/1e3)
    elif particles[i].label == 2: # ions
        x2.append(particles[i].x/1e-3)
        v2.append(particles[i].vx)

xe = []
for i in range(NBBins):
    xe.append(bins[i].x/1e-3)
E = []
for i in range(NBBins):
    E.append(bins[i].E)
phi = []
for i in range(NBBins):
    phi.append(bins[i].phi)
npos = []
for i in range(NBBins):
    npos.append(bins[i].npos*qsuper)
nneg = []
for i in range(NBBins):
    nneg.append(bins[i].nneg*qsuper)

# super particle positions
ax[0].set(ylabel='v_el (km/s), v_ion (m/s)',ylim=(-vscale/1e3,vscale/1e3),xlim=(xmin/1e-3,(xmax)/1e-3))
line1 = ax[0].plot(x1,v1,label='electrons',color='r',markersize=0.5,marker='o',lw=0)[0]
line2 = ax[0].plot(x2,v2,label='ions',color='b',markersize=0.5,marker='o',lw=0)[0]
ax[0].legend(loc=3,fontsize=4)

# info box
infobox = ''
infobox += '#particles: ' + "{:.0f}".format(NBParticles) + '\n'    
infobox += '#bins: ' + "{:.0f}".format(NBBins) + '\n'    
infobox += 'n0: ' + "{:.0e}".format(n0) + ' (1/m^3)\n'
infobox += 'Te: ' + "{:.1f}".format(Te) + ' (eV)\n'
infobox += 'l_Debye: ' + "{:.0e}".format(ldebye) + ' (m)\n'    
infobox += 'fpe: ' + "{:.0e}".format(fpe) + ' (1/s)\n'    
infobox += 'nu_el: ' + "{:.0e}".format(nuel) + ' (1/s)\n'    
infobox += 'nu_ion: ' + "{:.0e}".format(nuion) + ' (1/s)'    
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8) 
ax[0].text(0.05,0.95,infobox, fontsize=4,bbox=props,verticalalignment='top',transform=ax[0].transAxes)

# potential, Efield
axp = ax[0].twinx()
axp.set(ylabel='phi (V), E (kV/m)',ylim=(-phiscale-abs(V0),phiscale+abs(V0)))
linephi = axp.plot(xe,phi,label='phi',color='m',lw=2,alpha=0.8)[0]
lineE = axp.plot(xe,E,label='E',color='g',lw=2,alpha=0.8)[0]
axp.legend(loc=1,fontsize=4)

# densities
ax[1].set(xlabel='x (mm)',ylabel='n (m^-3)',ylim=(0,2*n0),xlim=(xmin/1e-3,(xmax)/1e-3))
linen1 = ax[1].plot(xe,nneg,label='electrons',color='r',lw=1)[0]
linen2 = ax[1].plot(xe,npos,label='ions',color='b',lw=1)[0]
ax[1].legend(loc=2,fontsize=4)

# energies for EEDF
enerbin = []
for i in range(NBEEDF):
    enerbin.append(EEDFmax/NBEEDF*i)
ener = []
for i in range(NBEEDF):
    ener.append(0)
ax[2].set(xlabel='energy (eV)',ylabel='f(E)',ylim=(1/NBParticles,0.5),xlim=(0,EEDFmax))
ax[2].set_yscale('log')
lineeedf = ax[2].plot(enerbin,ener,label='EEDF',color='g',lw=1)[0]

ERange = EEDFmax*1e-19
Etheory = np.linspace(1e-26,ERange,1000)
NEtheory = np.zeros(1000)
NEtheory2 = np.zeros(1000)
Ttheory = Te 
Ttheory2 = 3*Te
for i in range(1000):
        vt = np.sqrt(2*Etheory[i]/(me))
        NEtheory[i] = (me/(el*Ttheory))**(0.5)*np.exp((-0.5*me*vt**2)/(el*Ttheory))*1/(me*vt)*ERange/NBEEDF*1/np.pi
        NEtheory2[i] = (me/(el*Ttheory2))**(0.5)*np.exp((-0.5*me*vt**2)/(el*Ttheory2))*1/(me*vt)*ERange/NBEEDF*1/np.pi
lineetheory = ax[2].plot(Etheory/1e-19,NEtheory,color='r',linestyle='dashed',label='theory (Te:'+"{:.1f}".format(Ttheory)+')')[0]
lineetheory2 = ax[2].plot(Etheory/1e-19,NEtheory2,color='orange',linestyle='dashed',label='theory (Te:'+"{:.1f}".format(Ttheory2)+')')[0]
ax[2].legend(loc=1,fontsize=6)

t = 0
steps = int(tend/dtplot)


# -----------------------------------------------------------------------------
# Main Routine
# -----------------------------------------------------------------------------
def animate(k):
  global t,v,x,ax
  global xt,yt,zt,NBParticles
  
  time = 0
  while time < dtplot:    
    t += dt
    time += dt
    
    # 1st half kick 
    for i in range(NBParticles):
        particles[i].vx += particles[i].ax*0.5*dt
        particles[i].vy += particles[i].ay*0.5*dt

    # Update positions 
    for i in range(NBParticles):
        particles[i].x += particles[i].vx*dt
        # Particles lost at boundaries
        # then move last particle and replace it
        # for the i-th one
        if particles[i].x>=xmax or particles[i].x<=xmin:
            particles[i] = particles[NBParticles-1]
            particles[i].x += particles[i].vx*dt
            NBParticles -= 1   
        
    # PIC cycle 
    ParticlesToBins()
    CalculateEfield2(t)
    EFieldToParticles()
    ParticlesAcceleration(t)

    # ------------------------------------------------------
    # Collision calculation electron argon
    # from eduPIC code example by Z. Donko et al.
    # only electrons collide
    # ------------------------------------------------------
    for i in range(NBParticles):
      vx1 = particles[i].vx
      vy1 = particles[i].vy
      vz1 = particles[i].vz  
      if np.random.rand() < 1-np.exp(-nuel*dt) and particles[i].label == 1:
            # energy balance
            gx = vx1                             
            gy = vy1
            gz = vz1
            g  = np.sqrt(gx * gx + gy * gy + gz * gz)
            wx = massratio1 * vx1
            wy = massratio1 * vy1
            wz = massratio1 * vz1
            
            # calulate Euler angles
            if gx == 0:
                theta = 0.5*np.pi
            else:
                theta = np.arctan2(np.sqrt(gy**2+gz**2),gx)
            if gy == 0:
                if gz > 0:
                    phi = 0.5*np.pi
                else:
                    phi = -0.5*np.pi
            else:        
                phi = np.arctan2(gz,gy)

            st = np.sin(theta)
            ct = np.cos(theta)
            sp = np.sin(phi)
            cp = np.cos(phi)

            # select scattering angles 
            xi = np.arccos(1-2*np.random.rand())
            eta = 2*np.pi*np.random.rand()
        
            sc = np.sin(xi)
            cc = np.cos(xi)
            se = np.sin(eta)
            ce = np.cos(eta)
        
            # rotate velocity vector
            gx = g * (ct * cc - st * sc * ce)
            gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se)
            gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se)
      
            # velocity of the electron after collision in the lab system
            vx1 = wx + massratio2 * gx
            vy1 = wy + massratio2 * gy
            vz1 = wz + massratio2 * gz
    
      particles[i].vx = vx1 
      particles[i].vy = vy1
      particles[i].vz = vz1   

    # ------------------------------------------------------
    # Collision calculation argonion argon
    # from eduPIC code example by Z. Donko et al.
    # only electrons collide
    # ------------------------------------------------------
    for i in range(NBParticles):
      vx1 = particles[i].vx
      vy1 = particles[i].vy
      vz1 = particles[i].vz  
      if np.random.rand() < 1-np.exp(-nuion*dt) and particles[i].label == 2:
            # calculate relative velocity before collision
            # random Maxwellian target atom already selected 
            # (vx_2,vy_2,vz_2 velocity components of target atom come with the call)
            vx2 = 0
            vy2 = 0
            vz2 = 0
            gx = vx1-vx2;
            gy = vy1-vy2;
            gz = vz1-vz2;
            g  = np.sqrt(gx * gx + gy * gy + gz * gz)
            wx = 0.5 * (vx1 + vx2)
            wy = 0.5 * (vy1 + vy2)
            wz = 0.5 * (vz1 + vz2)

            # find Euler angles:
            if gx == 0:
                theta = 0.5 * np.pi 
            else:
                theta = np.arctan2(np.sqrt(gy * gy + gz * gz),gx)
            if (gy == 0):
                if (gz > 0):
                    phi = 0.5 * np.pi
                else:
                    phi = - 0.5 * np.pi
            else:
                phi = np.arctan2(gz, gy)


            # determine the type of collision based on cross sections and generate scattering angle

            chi = np.arccos(1.0 - 2.0 * np.random.rand())  # isotropic scattering angle
            eta = 2*np.pi* np.random.rand() #
            sc  = np.sin(chi)
            cc  = np.cos(chi)
            se  = np.sin(eta)
            ce  = np.cos(eta)
            st  = np.sin(theta)
            ct  = np.cos(theta)
            sp  = np.sin(phi)
            cp  = np.cos(phi)

            # compute new relative velocity:
            gx = g * (ct * cc - st * sc * ce);
            gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se)
            gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se)

            # post-collision velocity of the ion
            vx1 = wx + 0.5 * gx
            vy1 = wy + 0.5 * gy
            vz1 = wz + 0.5 * gz

      particles[i].vx = vx1 
      particles[i].vy = vy1
      particles[i].vz = vz1   
      

    # second half kick 
    for i in range(NBParticles):
        particles[i].vx += particles[i].ax*0.5*dt;
        particles[i].vy += particles[i].ay*0.5*dt;
     
    
  # ------------------------------------------------------
  # Display solution
  # ------------------------------------------------------    
  x1 = []
  v1 = []
  x2 = []
  v2 = []
  nelectrons = 0
  nions = 1
  for i in range(NBParticles):
      if particles[i].label == 1: # electrons
        x1.append(particles[i].x/1e-3)
        v1.append(particles[i].vx/1e3)
        nelectrons += 1
      elif particles[i].label == 2: # ions
        x2.append(particles[i].x/1e-3)
        v2.append(particles[i].vx)
        nions += 1
  E = []
  for i in range(NBBins):
    E.append(bins[i].E/1e3)
  phi = []
  for i in range(NBBins):
    phi.append(bins[i].phi)
  npos = []
  for i in range(NBBins):
    npos.append(bins[i].npos*qsuper)
  nneg = []
  for i in range(NBBins):
    nneg.append(bins[i].nneg*qsuper)
  
  ener = np.zeros(NBEEDF)  
  for i in range(NBParticles):
      if particles[i].label == 1:
         energy = particles[i].vx**2+particles[i].vy**2+particles[i].vz**2
         energy = energy*0.5*me/el
         #print(energy)
         enerbinID = int(energy/(EEDFmax/NBEEDF))
         #print(enerbinID)
         if enerbinID>=0 and enerbinID<NBEEDF:
             ener[enerbinID] += 1/NBParticles
  lineeedf.set_ydata(ener)
  
        
  line1.set_xdata(x1)
  line1.set_ydata(v1) 
  line2.set_xdata(x2)
  line2.set_ydata(v2)

  linephi.set_ydata(phi)
  lineE.set_ydata(E)
  
  linen1.set_ydata(nneg)
  linen2.set_ydata(npos)
        
  # Update the title of the plot
  fig.suptitle('time: ' + "{:.2f}".format(t/1e-9)+' ns')
  print('time: ' + "{:.2f}".format(t/1e-9)+' ns of '+ "{:.1f}".format(tend/1e-9)+' ns')
       
# ----------------------------------------
# Create animation
# ---------------------------------
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=180)

xp = []
vxp = []
vyp = []
qp = []
mp = []
labelp = []
for i in range(NBParticles):
    xp.append(particles[i].x)
    vxp.append(particles[i].vx)
    vyp.append(particles[i].vy)
    qp.append(particles[i].q)
    mp.append(particles[i].m)
    labelp.append(particles[i].label)
    
   
np.savetxt(initialfilename,[xp,vxp,vyp,qp,mp,labelp])    

np.savetxt('picrfoutel.txt', np.transpose([x1,v1]), fmt='%.3f', header="x1 v1")
np.savetxt('picrfoutE.txt', np.transpose([xe,E]), fmt='%.3f', header="x1 E")
np.savetxt('picrfoutphi.txt', np.transpose([xe,phi]), fmt='%.3f', header="x1 phi")
np.savetxt('picrfoutnpos.txt', np.transpose([xe,npos]), fmt='%.3f', header="x1 phi")
np.savetxt('picrfoutnneg.txt', np.transpose([xe,nneg]), fmt='%.3f', header="x1 phi")
                