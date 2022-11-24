# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 19:51:46 2022

@author: Achim
"""

import numpy as np
import matplotlib.pyplot as plt # for plot functons
import scipy as sp
import matplotlib.animation as animation


# --------------------------------------------------------------------------
# Define initial constants
# --------------------------------------------------------------------------
time = 0
t = 0
dtplot = 1e-9 # time distance in between frames for animation
simulationtime = 1000e-9 # total length of simulation in s
steps = int(simulationtime/dtplot) # total number of frames
dt=1e-12 # time constant for iteration


el = 1.6e-19
eps0 = 8.8e-12
me = 9e-31
kB = 1.38e-23
mp = 1.6e-27

mion = 40
p0 = 1 # pressure in Pa
tn = 300 # gas temperature
te0 = 3 # electron temperature in eV
zwidth = 0.1
V0 = -100

numdamp = 0 # numerical damping

znb = 500
z = np.zeros(znb+2)

ne = np.zeros(znb+2)
dne = np.zeros(znb+2)
ni = np.zeros(znb+2)
dni = np.zeros(znb+2)
Phi = np.zeros(znb+2)
EFeld = np.zeros(znb+2)
te = np.ones(znb+2)*te0
je = np.zeros(znb+2)
ji = np.zeros(znb+2)


a = np.zeros(znb+2)    
b = np.zeros(znb+2)    
c = np.zeros(znb+2)    
r = np.zeros(znb+2)    
ue = np.zeros(znb+2)    
ui = np.zeros(znb+2)    
   

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
def CalculatePotential():
   global a,b,c,r,ue,znb,EFeld,Phi 
   for i in range(1,znb+1):
      a[i]=1
      b[i]=-2
      c[i]=1
      r[i]=el/(eps0)*(ne[i]-ni[i])*(z[2]-z[1])**2
   nsurf = 0   
   r[1]=r[1]+el/(eps0)*nsurf*(z[2]-z[1])**2
   r[1]=r[1]+V0*(z[2]-z[1])**2
   MatInvTriagonal(a,b,c,r,ue,znb);
   for i in range(1,znb+1): 
      Phi[i]=ue[i]
   for i in range(1,znb+1):
      if i==1: 
         EFeld[i]=(Phi[i]-Phi[i+1])/(np.abs(z[2]-z[1]))
      elif i==znb: 
         EFeld[i]=(Phi[i-1]-Phi[i])/(np.abs(z[2]-z[1]))
      else:
         EFeld[i]=(Phi[i-1]-Phi[i+1])/(2*np.abs(z[2]-z[1]));

# --------------------------------------------------------------------------
# Initialize Distributions
# --------------------------------------------------------------------------
for i in range(1,znb+1):
    ne[i]=np.cos((i-1)/(znb-1)*np.pi-np.pi/2)*1e14+1e13;
    ni[i]=np.cos((i-1)/(znb-1)*np.pi-np.pi/2)*1e14+1e13;
    te[i]=te0;
    z[i]=-zwidth/2+(i-1)/(znb-1)*zwidth

#n_e_Zentrum=1e16 # m^-3
#Breite = 0.01 # m
#for i in range(1,znb+1):
#    te[i]=te0;
#    z[i]=-zwidth/2+(i-1)/(znb-1)*zwidth
#    ne[i] = n_e_Zentrum*np.exp(-z[i]**2/Breite**2)
#    ni[i] = n_e_Zentrum*np.exp(-z[i]**2/Breite**2)
    


ne[1]=0
ne[znb]=0
dz=(z[2]-z[1])


# --------------------------------------------------------------------------
# Define plots
# --------------------------------------------------------------------------
fig, ax = plt.subplots(2,1,figsize=(8,4.5))
line1 = ax[0].plot(z[1:znb],ne[1:znb],label='ne')[0]
line2 = ax[0].plot(z[1:znb],ni[1:znb],label='ni')[0] 
ax[0].set(xlabel='position',ylabel='density')
ax[0].legend()

line3 = ax[1].plot(z[1:znb],EFeld[1:znb],label='E field')[0]
ax[1].set(xlabel='position',ylabel='E field',ylim=(-5000,5000))
ax[1].legend() 

print('Start Simulation of '+"{:.0f}".format(steps)+' frames')

# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------
def animate(k):
  global t,ne,ni,dt,dne,dni,time

  print('time: '+ "{:.1f}".format(t/1e-9)+' ns of '+"{:.0f}".format(simulationtime/1e-9)+' ns' )
  while time<dtplot:
     time += dt   

     CalculatePotential()
     
     t += dt
     muel = mue(te[znb],tn,p0)
     Del = De(te[1],tn,p0)
     muion = mui(te[znb],tn,p0)
     Dion = Di(te[1],tn,p0)
     # -----------------------------------------------
     #  Inner Nodes
     #  Simple Lax-Wendroff scheme
     # ----------------------------------------------
     for i in range(2,znb):
         dne[i] = 0
         dni[i] = 0
     for i in range(2,znb):
        dne[i]=(-dt/(2*dz)*ne[i+1]*muel*(EFeld[i+1])
                +dt/(2*dz)*ne[i-1]*muel*(EFeld[i-1])
                +dt/(2*dz**2)*Del*(ne[i+1]-2*ne[i]+ne[i-1]))
        #dni[i]=(-dt/(2*dz)*ni[i+1]*muion*0.5*(EFeld[i+1]+EFeld[i])
        #        +dt/(2*dz)*ni[i-1]*muion*0.5*(EFeld[i-1]+EFeld[i])
        #        +dt/(2*dz**2)*Dion*(ni[i+1]-2*ni[i]+ni[i-1]))
        dni[i]=(-dt/(2*dz)*ni[i+1]*muion*(EFeld[i+1])
                +dt/(2*dz)*ni[i-1]*muion*(EFeld[i-1])
                +dt/(2*dz**2)*Dion*(ni[i+1]-2*ni[i]+ni[i-1]))
                
     #-------------------------------------------
     # BC electrons
     # ------------------------------------------}   
     dne[1]=-1/4*vth(te[1],tn,p0)*ne[1]*dt;
     dne[znb]=-1/4*vth(te[1],tn,p0)*ne[znb]*dt;
         
     # -------------------------------------------
     # BC ions
     #------------------------------------------}
     dni[1]=(1/((dz))*mui(te[1],tn,p0)*(EFeld[1])*ni[1]*dt
            -1/((dz))*mui(te[2],tn,p0)*(EFeld[2])*ni[2]*dt)
     dni[znb]=(-1/((dz))*mui(te[znb],tn,p0)*(EFeld[znb])*ni[znb]*dt
               +1/((dz))*mui(te[znb],tn,p0)*(EFeld[znb-1])*ni[znb-1]*dt)
     #dni[1] = 0 
     #dni[znb] = 0 
     
     # -----------------------------------------------
     # Propagate densities
     # ----------------------------------------------}
     for i in range(1,znb+1):
         ne[i] += dne[i];
         ni[i] += dni[i];
     #ni[znb] = ni[znb-1]    
     #ni[1] = ni[2]    
       
     # -----------------------------------------------
     # calculate fluxes
     # ----------------------------------------------}
     for i in range(1,znb+1):
        je[i] = mue(te[i],tn,p0)*ne[i]*EFeld[i];
        ji[i] = mui(te[i],tn,p0)*ni[i]*EFeld[i];
      
     line1.set_xdata(z[1:znb])
     line1.set_ydata(ne[1:znb]) 
     line2.set_xdata(z[1:znb])
     line2.set_ydata(ni[1:znb]) 
     line3.set_xdata(z[1:znb])
     line3.set_ydata(EFeld[1:znb]) 
     fig.suptitle('Time: ' + "{:.0f}".format(t/1e-9)+' ns')
  
  time = 0  

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save('driftdiffusionanimationexplicit.mp4',fps=25,dpi=180)

print('End Simulation')

#ax[0].plot(z[1:znb],ne[1:znb],label='ne')[0]
#ax[0].plot(z[1:znb],ni[1:znb],label='ni')[0] 
#ax[1].plot(z[1:znb],EFeld[1:znb],label='E field')[0]
