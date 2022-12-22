#
# Plasma oscillations
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
import matplotlib.pyplot as plt # for plot functons
import scipy as sp
import matplotlib.animation as animation


model = 'plasmafrequency'
model = 'debye'
#model = 'ionoscillation'
model = 'waveelectron'
model = 'waveion'


if model == 'plasmafrequency':
    dtplot = 1e-10 # time distance in between frames for animation
    simulationtime = 50e-9 # total length of simulation in s
    dt=1e-11 # time constant for iteration
    num =  1e8 # collision frequency
    animationfile = 'plasmafrequency_nu_8_v_1e4.mp4'
    n0 = 1e14 # charge density
    perturbation = 0.02 # perturbation of the quasineutrality
    dist = 'cosinuselectrons' # initial distribution
    efieldrange = 3e2 # range of efield plot in V/m
    vrange = 3e2 # range of velocity plot in km/s
    boundarycondition = 'default'
    zwidth = 0.1 # width of simulation volume in m
if model == 'debye':
    dtplot = 2e-11 # time distance in between frames for animation
    simulationtime = 10e-9 # total length of simulation in s
    dt=1e-11 # time constant for iteration
    num = 1e8  # 1e8 # collision frequency
    animationfile = 'debye.mp4'
    n0 = 1e14 # charge density
    perturbation = 0.02 # perturbation of the quasineutrality
    dist = 'gauss' # initial distribution
    efieldrange = 1e2 # range of efield plot in V/m
    vrange = 3e2 # range of velocity plot in km/s
    boundarycondition = 'default'
    zwidth = 0.1 # width of simulation volume in m
if model == 'ionoscillation':
    dtplot = 1e-9 # time distance in between frames for animation
    simulationtime = 1000e-9 #9 total length of simulation in s
    dt=1e-11 # time constant for iteration
    num = 0 # 1e8 # collision frequency
    animationfile = 'ionoscillation_nu_0.mp4'
    n0 = 1e14 # charge density
    perturbation = 0.02 # perturbation of the quasineutrality
    dist = 'cosinus' # initial distribution
    efieldrange = 3e2 # range of efield plot in V/m
    vrange = 3e2 # range of velocity plot in km/s
    boundarycondition = 'default'
    zwidth = 0.1 # width of simulation volume in m
if model == 'waveelectron':
    dtplot = 1e-10 # time distance in between frames for animation
    simulationtime = 100e-9 # total length of simulation in s
    dt=1e-11 # time constant for iteration
    num = 0  # 1e8 # collision frequency
    animationfile = 'waveelectron.mp4'
    n0 = 1e14 # charge density
    fexcitation = 5e7 #excitation v in central region
    aexcitation = 5
    perturbation = 0.02 # perturbation of the quasineutrality
    dist = 'wave' # initial distribution
    efieldrange = 1e2 # range of efield plot in V/m
    vrange = 3e2 # range of velocity plot in km/s
    boundarycondition = 'fixed'
    zwidth = 0.1 # width of simulation volume in m
if model == 'waveion':
    dtplot = 1e-9 # time distance in between frames for animation
    simulationtime = 100e-8 # total length of simulation in s
    dt=1e-11 # time constant for iteration
    num = 0  # 1e8 # collision frequency
    animationfile = 'waveion.mp4'
    fexcitation = 1e6
    aexcitation = 0.05
    n0 = 1e14 # charge density
    perturbation = 0.02 # perturbation of the quasineutrality
    dist = 'wave' # initial distribution
    efieldrange = 1e2 # range of efield plot in V/m
    vrange = 3e2 # range of velocity plot in km/s
    boundarycondition = 'fixed'
    zwidth = 0.1 # width of simulation volume in m
    
    
# --------------------------------------------------------------------------
# Define initial constants
# --------------------------------------------------------------------------
time = 0
t = 0
steps = int(simulationtime/dtplot) # total number of frames

el = 1.6e-19
eps0 = 8.8e-12
me = 9e-31
kB = 1.38e-23
mp = 1.6e-27

mion = 40*mp
p0 = 1 # pressure in Pa
tn = 300 # gas temperature
te0 = 3 # electron temperature in eV

V0 = 0

gammae = 1 # isothermal for electrons
Te = 3*el/kB # 

gammai = 5/3 # adiabatic for ions
Ti = 300


ldebye = np.sqrt(eps0*kB*Te/(n0*el**2))
wpe = np.sqrt(n0*el**2/(eps0*me))


znb = 500
z = np.zeros(znb+2)

ne = np.zeros(znb+2)
dne = np.zeros(znb+2)
ve = np.zeros(znb+2)
dve = np.zeros(znb+2)

ni = np.zeros(znb+2)
dni = np.zeros(znb+2)
vi = np.zeros(znb+2)
dvi = np.zeros(znb+2)


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
# Solving Poisson Equation 
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
# The initial distribution needs to be electrical neutral to
# allow the calculation of the electrical potential for periodic
# boundary condition
# --------------------------------------------------------------------------
if dist == 'cosinuselectrons':
    for i in range(1,znb+1):
        ne[i]=np.sin((i-1)/(znb)*8*np.pi)*perturbation*n0+n0
        ni[i]=n0;
        ve[i]=0
        vi[i]=0;
        te[i]=te0;
        z[i]=-zwidth/2+(i-1)/(znb)*zwidth
if dist == 'cosinus':
    for i in range(1,znb+1):
        ne[i]=np.sin((i-1)/(znb)*8*np.pi)*perturbation*n0+n0
        ni[i]=np.sin((i-1)/(znb)*8*np.pi)*perturbation*n0+n0
        ve[i]=0;
        vi[i]=0;
        te[i]=te0;
        z[i]=-zwidth/2+(i-1)/(znb)*zwidth
if dist =='wave':
    for i in range(0,znb+1):
        te[i]=te0;
        z[i]=-zwidth/2+(i-1)/(znb-1)*zwidth
        breitee = 0.002*zwidth
        breiteion = 0.02*zwidth
        # perturbation as time dependent gaussian trigger in velocity
        ne[i] = n0#+perturbation*n0*np.exp(-(z[i])**2/breitee**2)*1/(breitee)*breiteion
        ni[i] = n0      
        ve[i]=0;
        vi[i]=0;
elif dist =='gauss':
    for i in range(0,znb+1):
        te[i]=te0;
        z[i]=-zwidth/2+(i-1)/(znb-1)*zwidth
        breitee = 0.2*zwidth
        breiteion = 0.02*zwidth
        ne[i] = n0+perturbation*n0*np.exp(-(z[i])**2/breitee**2)*1/(breitee)*breiteion
        ni[i] = n0+perturbation*n0*np.exp(-(z[i])**2/breiteion**2)*1/(breiteion)*breiteion
        ve[i]=0;
        vi[i]=0;

# print(sum(ni)-sum(ne))    
dz=(z[2]-z[1])


# --------------------------------------------------------------------------
# Define plots
# --------------------------------------------------------------------------
fig, ax = plt.subplots(2,1,figsize=(8,4.5))
plt.subplots_adjust(right=0.85)

ax[0].set(xlabel='position (m)',ylabel='density (1/m^3)',ylim=(n0-(perturbation+0.01)*n0,n0+(perturbation+0.01)*n0),xlim=(-zwidth/2,zwidth/2))
linene = ax[0].plot(z[1:znb],ne[1:znb],label='ne',color='r')[0]
lineni = ax[0].plot(z[1:znb],ni[1:znb],label='ni',color='b')[0] 
ax[0].legend(loc=4,fontsize=6)

axp1 = ax[0].twinx()
axp1.set(xlabel='position (m)',ylabel='velocity (km/s)',ylim=(-vrange,vrange),label='ve')
lineve = axp1.plot(z[1:znb],ve[1:znb]/1e3,label='ve',linestyle='dashed',color='r')[0]
linevi = axp1.plot(z[1:znb],vi[1:znb]/1e3,label='vi',linestyle='dashed',color='b')[0] 
axp1.legend(loc=1,fontsize=6)

linee = ax[1].plot(z[1:znb],EFeld[1:znb],label='E field')[0]
ax[1].set(xlabel='position (m)',ylabel='E field (V/m)',ylim=(-efieldrange,efieldrange),xlim=(-zwidth/2,zwidth/2))
ax[1].legend(loc=2,fontsize=6) 

# info box
infobox = ''
infobox += 'nu: ' + "{:.0e}".format(num) + ' (1/s) \n'   
infobox += 'l_Debye: ' + "{:.0e}".format(ldebye) + ' (m) \n'   
infobox += 'w_plasma,e: ' + "{:.0e}".format(wpe) + ' (1/s)'   
#if model == "electronoscillation":
#    infobox += 'nu: ' + "{:.0e}".format(num) + ' (1/s)'   
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[0].text(0.02,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax[0].transAxes)


print('Start Simulation of '+"{:.0f}".format(steps)+' frames')

# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------
def animate(k):
  global t,ne,ni,dt,dne,dni,time

  print('time: '+ "{:.2f}".format(t/1e-9)+' ns of '+"{:.0f}".format(simulationtime/1e-9)+' ns ')
  while time<dtplot:
     time += dt   

     CalculatePotential()
     
     t += dt
     
     #if dist=='wave':
     #    zpt = int(znb/2)
     #    ne[zpt] = n0+perturbation*n0*np.sin(t*1e8)
         
     # -----------------------------------------------
     #  Inner Nodes
     # ----------------------------------------------
     for i in range(2,znb):
        q = -1
        neleft = 0.5*(ne[i-1]+ne[i])
        neright = 0.5*(ne[i]+ne[i+1]) 
        veleft = 0.5*(ve[i-1]+ve[i])
        veright = 0.5*(ve[i]+ve[i+1]) 
        Eleft = 0.5*(EFeld[i-1]+EFeld[i])
        Eright = 0.5*(EFeld[i+1]+EFeld[i])
        #dne[i]=(-dt/(dz)*(neright*veright-neleft*veleft))
        dne[i]=(-dt/(dz)*(ne[i]*(veright-veleft)+ve[i]*(neright-neleft)))
        dve[i]=(dt*q*el/me*(EFeld[i]) - gammae*kB*Te/(me*ne[i])*(neright-neleft)/(dz)*dt
                 -ve[i]*num*dt-ve[i]*(veright-veleft)/dz*dt)
        
        q = 1
        nileft = 0.5*(ni[i-1]+ni[i])
        niright = 0.5*(ni[i]+ni[i+1]) 
        vileft = 0.5*(vi[i-1]+vi[i])
        viright = 0.5*(vi[i]+vi[i+1]) 
        #dni[i]=(-dt/(dz)*(niright*viright-nileft*vileft))
        dni[i]=(-dt/(dz)*(ni[i]*(viright-vileft)+vi[i]*(niright-nileft)))
        dvi[i]=(dt*q*el/mion*(EFeld[i]) - gammai*kB*Ti/(mion*ni[i])*(niright-nileft)/(dz)*dt
                 -vi[i]*num*dt-vi[i]*(viright-vileft)/dz*dt)


        # add perturbation in ve or vi for wave models        
        if model == 'waveelectron':
          dve[i] += np.sin(t*fexcitation*2*np.pi)*aexcitation*np.exp(-(z[i])**2/breitee**2)*1/(breitee)
        if model == 'waveion':
          dvi[i] += np.sin(t*fexcitation*2*np.pi)*aexcitation*np.exp(-(z[i])**2/breiteion**2)*1/(breiteion)
        
     if boundarycondition == 'fixed':   
       #-------------------------------------------
       # BC left
       # ------------------------------------------}   
       dne[1]=0
       dve[1]=dve[2]
       dni[1]=0
       dvi[1]=dvi[2]
     
       #-------------------------------------------
       # BC right
       # ------------------------------------------}   
       dne[znb]=0
       dve[znb]=dve[znb-1]
       dni[znb]=0
       dvi[znb]=dvi[znb-1]             
     else:   
       #-------------------------------------------
       # BC left
       # ------------------------------------------}   
       #BC electrons
       q = -1
       neleft = 0.5*(ne[znb]+ne[1])
       neright = 0.5*(ne[1]+ne[2]) 
       veleft = 0.5*(ve[znb]+ve[1])
       veright = 0.5*(ve[1]+ve[2]) 
       dne[1]=(-dt/(dz)*(neleft*veleft-neright*veright)*(-1))
       dve[1]=(dt*q*el/me*EFeld[1] - gammae*kB*Te/(me*ne[1])*(neright-neleft)/(dz)*dt
                 -ve[1]*num*dt-ve[1]*(veleft-veright)*(-1)/dz*dt)
       #BC ions
       q = 1
       nileft = 0.5*(ni[znb]+ni[1])
       niright = 0.5*(ni[1]+ni[2]) 
       vileft = 0.5*(vi[znb]+vi[1])
       viright = 0.5*(vi[1]+vi[2]) 
       dni[1]=(-dt/(dz)*(nileft*vileft-niright*viright)*(-1))
       dvi[1]=(dt*q*el/mion*EFeld[1] - gammai*kB*Ti/(mion*ni[1])*(niright-nileft)/(dz)*dt
                 -vi[1]*num*dt-ve[1]*(vileft-viright)*(-1)/dz*dt)

     
       #-------------------------------------------
       # BC right
       # ------------------------------------------}   
       #BC electrons
       q = -1
       neleft = 0.5*(ne[znb-1]+ne[znb])
       neright = 0.5*(ne[znb]+ne[1]) 
       veleft = 0.5*(ve[znb-1]+ve[znb])
       veright = 0.5*(ve[znb]+ve[1]) 
       dne[znb]=(-dt/(dz)*(neleft*veleft-neright*veright)*(-1))
       dve[znb]=(dt*q*el/me*EFeld[znb] - gammae*kB*Te/(me*ne[znb])*(neright-neleft)/(dz)*dt
                 -ve[znb]*num*dt-ve[znb]*(veleft-veright)*(-1)/dz*dt)
       #BC ions
       q = 1
       nileft = 0.5*(ni[znb-1]+ni[znb])
       niright = 0.5*(ni[znb]+ni[1]) 
       vileft = 0.5*(vi[znb-1]+vi[znb])
       viright = 0.5*(vi[znb]+vi[1]) 
       dni[znb]=(-dt/(dz)*(nileft*vileft-niright*viright)*(-1))
       dvi[znb]=(dt*q*el/mion*EFeld[znb] - gammai*kB*Ti/(mion*ni[znb])*(niright-nileft)/(dz)*dt
                 -vi[znb]*num*dt-vi[znb]*(vileft-viright)*(-1)/dz*dt)

     
     # -----------------------------------------------
     # Propagate densities
     # ----------------------------------------------}
     for i in range(1,znb+1):
         #ve[i] += np.sin(t*1e7)*1000*np.exp(-(z[i])**2/breitee**2)*1/(breitee)*breiteion
         ne[i] += dne[i];
         ni[i] += dni[i];
         ve[i] += dve[i];
         vi[i] += dvi[i];
       
      
     linene.set_xdata(z[1:znb])
     linene.set_ydata(ne[1:znb]) 
     lineve.set_xdata(z[1:znb])
     lineve.set_ydata(ve[1:znb]/1e3) 
 
     lineni.set_xdata(z[1:znb])
     lineni.set_ydata(ni[1:znb]) 
     linevi.set_xdata(z[1:znb])
     linevi.set_ydata(vi[1:znb]/1e3) 

     linee.set_xdata(z[1:znb])
     linee.set_ydata(EFeld[1:znb]) 
     fig.suptitle('Time: ' + "{:.0f}".format(t/1e-9)+' ns')
  
  time = 0 
  #print(sum(ne-ni))

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationfile,fps=25,dpi=300)

print('End Simulation')

