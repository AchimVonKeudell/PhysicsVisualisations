#
# Drift diffusion approximation
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.animation as animation


# ----------------------------------------------------------------------------
# Setting up models
# ----------------------------------------------------------------------------

#model = 'RFCycle' # several RF cycles
model = 'ProfileNegativeIons' # Profileevolution DC
#model = 'ProfileEvolutionDC' # Profileevolution DC voltage
#model = 'ProfileGauss' # Profileevolution DC
#model = 'ProfileDecay' # Profile Decay
# model = 'Test' # Profile Decay

if model == 'RFCycle':
    V0 = 100
    dist = 'cosinus' # initial distribution
    dtplot = 1e-9
    simulationtime = 20*74e-9
    dt=1e-10
    simulationfile = 'RFCycleAnimation.gif'
if model == 'ProfileNegativeIons':
    V0 = 0
    dist = 'cosinus' # initial distribution
    dtplot = 10e-8
    simulationtime = 20000e-9
    alpha = 0.5 # fraction negative ions
    dt=2e-10
    simulationfile = 'ProfileNegativeIons.gif'
if model == 'ProfileEvolutionDC':
    V0 = 100
    dist = 'const' # initial distribution
    dtplot = 1e-10
    simulationtime = 200e-9
    dt=0.5e-10
    simulationfile = 'ProfileDCAnimation.gif'
if model == 'ProfileGauss':
    V0 = 0
    dist = 'gauss' # initial distribution
    dtplot = 2e-10
    simulationtime = 20e-9
    dt=1e-10
    simulationfile = 'ProfileGaussAnimation.gif'
if model == 'ProfileDecay':
    V0 = 0
    dist = 'decay' # initial distribution
    dtplot = 2e-9
    simulationtime = 10000e-9
    dt=2e-10
    simulationfile = 'ProfileDecayAnimation.mp4'
if model == 'Test':
    V0 = 100
    dist = 'decay' # initial distribution
    dtplot = 2e-9
    simulationtime = 5e-9
    dt=2e-10
    simulationfile = 'ProfileTest.mp4'


# ----------------------------------------------------------------------------
# General parameter
# ----------------------------------------------------------------------------

steps = int(simulationtime/dtplot)

el = 1.6e-19
eps0 = 8.8e-12
me = 9e-31
kB = 1.38e-23
mp = 1.6e-27

mion = 40
p0 = 1 # pressure in Pa
tn = 300 # gas temperature
te0 = 3 # electron temperature in eV
fRF = 13.56e6

numdamp = 0.005 # numerical damping

znb = 800
ne = np.zeros(znb+2)
ni = np.zeros(znb+2)
nn = np.zeros(znb+2)
Phi = np.zeros(znb+2)
EFeld = np.zeros(znb+2)
te = np.ones(znb+2)*te0
je = np.zeros(znb+2)
ji = np.zeros(znb+2)

z = np.zeros(znb+2)
zwidth = 0.1
dz=zwidth/znb

a = np.zeros(znb+2)    
b = np.zeros(znb+2)    
c = np.zeros(znb+2)    
r = np.zeros(znb+2)    
ue = np.zeros(znb+2)    
ui = np.zeros(znb+2)    
un = np.zeros(znb+2)    
   
# ----------------------------------------------------------------------------
# Voltage
# ----------------------------------------------------------------------------
def voltage(tt):
    if model=='RFCycle': # several RF cycles
        return V0*1*(np.cos(2*np.pi*fRF*tt))
    if model == 'ProfileNegativeIons': # Profileevolution DC
        return 0
    if model == 'ProfileEvolutionDC': # Profileevolution DC voltage
        return V0
    if model == 'ProfileDecay': # Profile Decay
        return V0
    
    
# ----------------------------------------------------------------------------
# Inverse of a Trigonal Matrix from Numerical Recipes
#
# | b1  c1  0   ... | u1   r1
# | a2  b2  c2  ... | u2 = r2
# | 0   a3  b3   c3 | u3 = r3
# | 0   0   a4   b4 | u4 = r4
# Solves for vector u
# ----------------------------------------------------------------------------   
def MatInvTriagonal(a,b,c,r,u,n):
   bet = b[1]
   u[1] = r[1]/bet
   gam = np.zeros(n+1)
   for j in range(2,n): 
      gam[j] = c[j-1]/bet
      bet = b[j]-a[j]*gam[j]
      u[j] = (r[j]-a[j]*u[j-1])/bet
   for j in range(n-1,0,-1): 
      u[j] = u[j]-gam[j+1]*u[j+1]
   return a,b,c,r,u


# ----------------------------------------------------------------------------
# Physical Parameter
# ----------------------------------------------------------------------------   
# Stossfrequenzen
def vcollisionelectron(te,tn,p):
  return (p/(kB*tn))*5.3e-13

def vcollisionion(te,tn,p):
  return (p/(kB*tn))*8.3e-16
 
# Diffusionskonstante der Elektronen
# De := k Te / (m vcollision)
def De(te,tn,p):
  return (te*el)/(me*vcollisionelectron(te,tn,p))

# thermische Geschw. der Elektronen
# vth = sqrt(e/(m vcollsion)
def vth(te,tn,p):
  return np.sqrt(3*el*te/(me))

# Beweglichkeit der Elektronen
# mue = e/(m vcollsion)
def mue(te,tn,p):
  return -el/(me*vcollisionelectron(te,tn,p));

# Diffusionskonstante der Ionen
# Di := k Ti / (m vcollision)
def Di(te,tn,p):
  return (kB*tn)/(mion*mp*vcollisionion(te,tn,p))

# Beweglichkeit der Ionen
# mui = e/(m vcollsion)
def mui(te,tn,p):
  return el/(mion*mp*vcollisionion(te,tn,p));

# Ionizationsfrequenz
def vionization(te,tn,p):
  return 0  
  # return 4e-21*np.sqrt(8*te/(np.pi*9.1e-31))*(p/(1.23e-23*tn))*np.exp(-24.6/te);
   
    
    
# ----------------------------------------------------------------------------   
# Solving Poissons Equation 
#----------------------------------------------------------------------------
def CalculatePotential(tvoltage):
   global a,b,c,r,ue,znb,EFeld,Phi,dz 
   for i in range(1,znb+1):
      a[i]=1
      b[i]=-2
      c[i]=1
      r[i]=el/(eps0)*(ne[i]-ni[i]+nn[i])*(dz)**2 
   #nsurf = 0    
   #r[1]=r[1]+el/(eps0)*nsurf*(z[2]-z[1])**2
   r[1]=voltage(tvoltage)
   r[znb] = 0
   
   MatInvTriagonal(a,b,c,r,ue,znb);
   for i in range(1,znb+1): 
      Phi[i]=ue[i]
   for i in range(2,znb+1):
      EFeld[i]=(Phi[i-1]-Phi[i+1])/(2*dz);
   EFeld[1]=(Phi[1]-Phi[2])/(dz)
   EFeld[znb]=(Phi[znb-1]-Phi[znb])/(dz)
 

# --------------------------------------------------------------------------
# Initialize Distributions Cosinus
# --------------------------------------------------------------------------
if dist == 'cosinus':
  for i in range(1,znb+2):
    ne[i]=(1-alpha)*(np.cos((i-1)/(znb-1)*np.pi-np.pi/2)*1e14+1e13)
    ni[i]=1.0*(np.cos((i-1)/(znb-1)*np.pi-np.pi/2)*1e14+1e13)
    nn[i]=alpha*(np.cos((i-1)/(znb-1)*np.pi-np.pi/2)*1e14+1e13)
    te[i]=te0;
    z[i]=-zwidth/2+(i-1)*dz
if dist == 'const':
  for i in range(1,znb+2):
    ne[i]=0.9*(1e14+1e13)
    ni[i]=1e14+1e13;
    nn[i]=0.1*(1e14+1e13)
    te[i]=te0;
    z[i]=-zwidth/2+(i-1)*dz
if dist == 'gauss':    
  n_e_Zentrum=1e14 # m^-3
  Breite = 0.01 # m
  for i in range(1,znb+1):
     te[i]=te0;
     z[i]=-zwidth/2+(i-1)/(znb-1)*zwidth
     ne[i] = 0.9*(n_e_Zentrum*np.exp(-z[i]**2/Breite**2))
     ni[i] = n_e_Zentrum*np.exp(-z[i]**2/Breite**2)
     nn[i] = 0.1*(n_e_Zentrum*np.exp(-z[i]**2/Breite**2))
if dist == 'decay':
  for i in range(1,znb+2):
    ne[i]=0.9*(np.cos((i-1)/(znb-1)*2*np.pi-np.pi/2)**2*1e14+1e13)
    ni[i]=np.cos((i-1)/(znb-1)*2*np.pi-np.pi/2)**2*1e14+1e13
    nn[i]=0.1*(np.cos((i-1)/(znb-1)*2*np.pi-np.pi/2)**2*1e14+1e13)    
    te[i]=te0;
    z[i]=-zwidth/2+(i-1)*dz
    


ne[1]=0
ne[znb]=0

t = 0
time = 0

lb = 2
ub = znb-1

fig, ax = plt.subplots(2,1,figsize=(8,4.5))
plt.subplots_adjust(right=0.85)

line1 = ax[0].plot(z[lb:ub],ne[lb:ub],label='ne',color='r')[0]
line2 = ax[0].plot(z[lb:ub],ni[lb:ub],label='ni',color='b')[0] 
line2n = ax[0].plot(z[lb:ub],nn[lb:ub],label='nn',color='g')[0] 
ax[0].set(xlabel='position (m)',ylabel='density (m^-3)',xlim=(-zwidth/2,zwidth/2),ylim=(0,1.3e14))
ax[0].legend(loc=0)

line3 = ax[1].plot(z[lb:ub],EFeld[lb:ub]/1e3,label='E field',color='g')[0]
ax[1].set(xlabel='position (m)',ylabel='E field (kV/m)',xlim=(-zwidth/2,zwidth/2),ylim=(-10,10))
ax[1].legend(loc=2)

ax2 = ax[1].twinx()
line4 = ax2.plot(z[lb:ub],Phi[lb:ub],label='Potential',color='m')[0]
if V0 != 0:
    ax2.set(ylim=(-1.5*V0,1.5*V0))
else:
    ax2.set(ylim=(-100,100))
    
ax2.set_ylabel('Potential (V)') 
ax2.legend(loc=4)

print('Start Simulation of '+"{:.0f}".format(steps)+' frames')


# --------------------------------------------------------------------------
#
# Implict Scheme using Lax method for treating the convection term
#
# Explicit formula (in the implicit case all values on the rhs ar at t+dt)
#
# (n_i^{t+dt} - n_i^t)/dt =
#
# Dleft (n_{i-1} - n_i)/dx^2 + Dleft (n_{i-1} - n_i)/dx^2
#
# -mue Eright n_{i+1}/(2 dx) + mue Eleft n_{i-1}/(2 dx)
#
# -------------------------------------------------------------------------- 

def animate(k):
  global t,ne,ni,dt,time

  print('time: '+ "{:.1f}".format(t/1e-9)+' ns of '+"{:.0f}".format(simulationtime/1e-9)+' ns' )
    
  while time<dtplot:     
     time += dt
     CalculatePotential(t)
     t += dt
     # -----------------------------------------------
     #  electrons
     # ----------------------------------------------
     for i in range(2,znb):
        p1=1/dt
        
        p21=1/(2*dz)*mue(te[i],tn,p0)*(EFeld[i-1])
        p22=1/(2*dz)*mue(te[i],tn,p0)*(EFeld[i+1])
        
        p3=1/dz**2*De(te[i],tn,p0)
        p5=vionization(te[i],tn,p0)

        a[i]=(-p3-p21)
        b[i]=(p1+2*p3)
        c[i]=(+p22-p3)
        r[i]=p1*ne[i]+p5*ne[i]

     #-------------------------------------------
     # Randbedingung Elektronen fuer z = 1
     # ------------------------------------------}
     p1=1/dt
    
     p21=1/(2*dz)*mue(te[2],tn,p0)*(EFeld[2])
     p3=1/dz**2*De(te[1],tn,p0)
    
     plb=1/4*vth(te[1],tn,p0)
     p5=vionization(te[1],tn,p0)

     a[1]=0
     b[1]=p1-plb+2*p3
     c[1]=-p3+p21
     r[1]=p1*ne[1]+p5*ne[1]


     # -----------------------------------------
     # Randbedingung Elektronen fuer z = znb
     # ------------------------------------------}
     p1=1/dt;
   
     p22=1/(2*dz)*mue(te[znb-1],tn,p0)*(EFeld[znb-1])
     p3=1/dz**2*De(te[znb],tn,p0);
    
     prb=1/(4)*vth(te[znb],tn,p0);
     p5=vionization(te[znb],tn,p0);

     a[znb]=-p3-p21;
     b[znb]=p1-prb+2*p3;
     c[znb]=0;
     r[znb]=p1*ne[znb]+p5*ne[znb];


     # -----------------------------------------
     # Invert matrix
     # ------------------------------------------}
  
     MatInvTriagonal(a,b,c,r,ue,znb);

     # -----------------------------------------------
     #  Neue Ionen-Dichten ausrechnen
     #--------------------------------------------}
     for i in range(2,znb): 
        p1=1/dt;
  
        p21=1/(2*dz)*mui(te[i],tn,p0)*(EFeld[i-1]);
        p22=1/(2*dz)*mui(te[i],tn,p0)*(EFeld[i+1]);
       
        p3=1/(dz**2)*Di(te[i],tn,p0);
        p5=vionization(te[i],tn,p0);
        
        vel0=np.abs(mue(te[i],tn,p0)*EFeld[i]);
        p4=numdamp*vel0/(dz);

        a[i]=(-p3-p21-p4);
        b[i]=(p1+2*p3+2*p4);
        c[i]=(p22-p3-p4);
        r[i]=p1*ni[i]+p5*ne[i];

     # -------------------------------------------
     # BC ions left edge
     #------------------------------------------}
     p1=1/dt;
          
     p21=1/((2*dz))*mui(te[2],tn,p0)*(EFeld[2])
     plb=1/((2*dz))*mui(te[1],tn,p0)*EFeld[1];      
     p3=1/(dz**2)*Di(te[1],tn,p0);
             
     p5=vionization(te[1],tn,p0);

     a[1]=0;
     b[1]=p1+2*p3-plb;
     c[1]=-p3+p21;
     r[1]=p1*ni[1]+p5*ne[1];


     #-----------------------------------------
     # BC ions right edge
     #-----------------------------------------}
     p1=1/dt;

     p21=1/((2*dz))*mui(te[znb-1],tn,p0)*(EFeld[znb-1])
     prb=1/((2*dz))*mui(te[znb],tn,p0)*EFeld[znb];      
     p3=1/(dz**2)*Di(te[znb],tn,p0);
     
     
     p5=vionization(te[znb],tn,p0);
     
     vel0=np.abs(mue(te[i],tn,p0)*EFeld[znb]);
     p4=numdamp*vel0/(dz);

     a[znb-1]=-p3-p21-p4;
     b[znb]=p1+2*p3-prb+2*p4;
     c[znb]=0;
     r[znb]=p1*ni[znb]+p5*ne[znb];

     
     # -----------------------------------------
     # Invert matrix
     # ------------------------------------------} 
     MatInvTriagonal(a,b,c,r,ui,znb);

    
     # -----------------------------------------------
     #  Neue negative Ionen-Dichten ausrechnen
     #--------------------------------------------}
     for i in range(2,znb): 
        p1=1/dt;
  
        p21=1/(2*dz)*(-1)*mui(te[i],tn,p0)*(EFeld[i-1]);
        p22=1/(2*dz)*(-1)*mui(te[i],tn,p0)*(EFeld[i+1]);
       
        p3=1/(dz**2)*Di(te[i],tn,p0);
        p5=vionization(te[i],tn,p0);
        
        vel0=np.abs(mue(te[i],tn,p0)*EFeld[i]);
        p4=numdamp*vel0/(dz);

        a[i]=(-p3-p21-p4);
        b[i]=(p1+2*p3+2*p4);
        c[i]=(p22-p3-p4);
        r[i]=p1*nn[i]+p5*ne[i];

     # -------------------------------------------
     # BC neg ions left edge
     #------------------------------------------}
     p1=1/dt;
          
     p21=1/((2*dz))*(-1)*mui(te[2],tn,p0)*(EFeld[2])
     plb=1/((2*dz))*(-1)*mui(te[1],tn,p0)*EFeld[1];      
     p3=1/(dz**2)*Di(te[1],tn,p0);
             
     p5=vionization(te[1],tn,p0);

     a[1]=0;
     b[1]=p1+2*p3-plb;
     c[1]=-p3+p21;
     r[1]=p1*nn[1]+p5*ne[1];


     #-----------------------------------------
     # BC ions right edge
     #-----------------------------------------}
     p1=1/dt;

     p21=1/((2*dz))*(-1)*mui(te[znb-1],tn,p0)*(EFeld[znb-1])
     prb=1/((2*dz))*(-1)*mui(te[znb],tn,p0)*EFeld[znb];      
     p3=1/(dz**2)*Di(te[znb],tn,p0);
     
     
     p5=vionization(te[znb],tn,p0);
     
     vel0=np.abs(mue(te[i],tn,p0)*EFeld[znb]);
     p4=numdamp*vel0/(dz);

     a[znb-1]=-p3-p21-p4;
     b[znb]=p1+2*p3-prb+2*p4;
     c[znb]=0;
     r[znb]=p1*nn[znb]+p5*ne[znb];

     
     # -----------------------------------------
     # Invert matrix
     # ------------------------------------------} 
     MatInvTriagonal(a,b,c,r,un,znb);

     # -----------------------------------------------
     # Neue Dichten uebergeben
     # ----------------------------------------------}
     sum=0;
     for i in range(1,znb+2):
         ne[i] = ue[i];
         sum += ne[i];
         ni[i] = ui[i];    
         nn[i] = un[i];    
    
    
     # -----------------------------------------------
     # calculate fluxes
     # ----------------------------------------------}
     for i in range(1,znb):
        je[i] = mue(te[i],tn,p0)*ne[i]*EFeld[i];
        ji[i] = mui(te[i],tn,p0)*ni[i]*EFeld[i];
      
     line1.set_xdata(z[lb:ub])
     line1.set_ydata(ne[lb:ub]) 
     line2.set_xdata(z[lb:ub])
     line2.set_ydata(ni[lb:ub]) 
     line2n.set_xdata(z[lb:ub])
     line2n.set_ydata(nn[lb:ub]) 
     
     line3.set_xdata(z[lb:ub])
     line3.set_ydata(EFeld[lb:ub]/1e3) 
     line4.set_xdata(z[lb:ub])
     line4.set_ydata(Phi[lb:ub]) 
     fig.suptitle('Time: ' + "{:.0f}".format(t/1e-9)+' ns')
  
  time = 0  


anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(simulationfile,fps=25,dpi=180)

#ax[0].plot(z[1:znb],ne[1:znb],color='r')
#ax[0].plot(z[1:znb],ni[1:znb],color='b')
#ax[1].plot(z[1:znb],EFeld[1:znb]/1e3)
#ax2.plot(z[1:znb],Phi[1:znb])
