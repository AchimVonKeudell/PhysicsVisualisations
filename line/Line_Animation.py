#
# Transmission line
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
import matplotlib.pyplot as plt # for plot functons
import scipy as sp
import matplotlib.animation as animation


model = 'fixed boundaries'
model = 'absorbing boundaries'
#model = 'open boundaries'

dtplot = 2e-10 # time distance in between frames for animation
simulationtime = 100e-9 # total length of simulation in s
dt=1e-12 # time constant for iteration
RefractiveIndexMedium = 1.5
RefractiveIndexAmbient = 1
RefractiveIndexPlug = 3
boundary1 = 11
boundary2 = 12.2
gamma = 0 # damping in the medium
    
if model == 'fixed boundaries':
    animationfile = 'transmission_fixed.mp4'
    BoundaryConditionLeft = 'fixed'
    BoundaryConditionRight = 'fixed'
if model == 'open boundaries':
    BoundaryConditionLeft = 'open'
    BoundaryConditionRight = 'open'
    animationfile = 'transmission_open.mp4'
if model == 'absorbing boundaries':
    BoundaryConditionLeft = 'absorbing'
    BoundaryConditionRight = 'absorbing'
    animationfile = 'transmission_absorbing.mp4'
    
# --------------------------------------------------------------------------
# Define initial constants
# --------------------------------------------------------------------------
time = 0
t = 0
steps = int(simulationtime/dtplot) # total number of frames

Xsize = 10000
dx = 0.002

line = np.zeros([Xsize+1,3])
x = np.linspace(0,dx*Xsize,Xsize+1)
dline = np.zeros(Xsize+1)
alphavector = np.zeros(Xsize+1)
signal = np.zeros(Xsize+1)
signalreflected = np.zeros(Xsize+1)
  

# --------------------------------------------------------------------------
# Initialize Distributions
# The initial distribution needs to be electrical neutral to
# allow the calculation of the electrical potential for periodic
# boundary condition
# --------------------------------------------------------------------------

# Starting Condition Amplitude 
for i in range(Xsize): 
    line[i,0] = np.exp(-(i-Xsize/4)**2/(0.3/dx)**2)
   
vphase = 3e8
gamma0 = 0
alpha = vphase*dt/dx

if alpha>1: print('von Neumann condition violated')
  

for i in range(Xsize):
      if i<=boundary1/dx: 
        alphavector[i] = alpha
      else:
        alphavector[i] = alpha*1/(RefractiveIndexMedium)
   
    
# Starting Condition Velocity 
for i in range(2,Xsize-1):
      dline[i] = -vphase/RefractiveIndexAmbient*(line[i+1,0]-line[i-1,0])/(2*dx)
    

# Starting Condition Ghost Point in Time }
for i in range(Xsize-1):
      line[i,1] = line[i,0]+0.5*alpha**2*(line[i+1,0]-2*line[i,0]+line[i-1,0]) + dline[i]*dt;
  
 
# --------------------------------------------------------------------------
# Define plots
# --------------------------------------------------------------------------
fig, ax = plt.subplots(1,1,figsize=(8,4.5))
plt.subplots_adjust(right=0.85)
yrange = 2

ax.set(xlabel='position (m)',ylabel='amplitude',ylim=(-yrange,yrange),xlim=(0,Xsize*dx))
lineamp = ax.plot(x,line[:,0],color='b')[0]
b3x = [boundary1,boundary2]
b3y = [yrange,yrange]
b3y2 = [-yrange,-yrange]
ax.fill_between(b3x,b3y,b3y2,color='lightblue',lw=1,alpha=0.8) 
b4x = [boundary2,Xsize*dx]
b4y = [yrange,yrange]
b4y2 = [-yrange,-yrange]
ax.fill_between(b4x,b4y,b4y2,color='orange',lw=1,alpha=0.8) 


ax.legend(loc=3,fontsize=6)

#axp1 = ax[0].twinx()
#axp1.set(xlabel='position (m)',ylabel='velocity (km/s)',ylim=(-vrange,vrange),label='ve')
#lineve = axp1.plot(z[1:znb],ve[1:znb]/1e3,label='ve',linestyle='dashed',color='r')[0]
#linevi = axp1.plot(z[1:znb],vi[1:znb]/1e3,label='vi',linestyle='dashed',color='b')[0] 
#axp1.legend(loc=1,fontsize=6)

#linee = ax[1].plot(z[1:znb],EFeld[1:znb],label='E field')[0]
#ax[1].set(xlabel='position (m)',ylabel='E field (V/m)',ylim=(-efieldrange,efieldrange),xlim=(-zwidth/2,zwidth/2))
#ax[1].legend(loc=2,fontsize=6) 

# info box
infobox = ''
infobox += 'model: ' + model
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)


infobox = ''
infobox += 'n ambient: ' + "{:.1f}".format(RefractiveIndexAmbient)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.35,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)

infobox = ''
infobox += 'n interface: ' + "{:.1f}".format(RefractiveIndexPlug)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(boundary1/(Xsize*dx)+0.02,0.75,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)

infobox = ''
infobox += 'n medium: ' + "{:.1f}".format(RefractiveIndexMedium)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(boundary2/(Xsize*dx)+0.02,0.35,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)


print('Start Simulation of '+"{:.1f}".format(steps)+' frames')

# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------
t = 0
def animate(k):
  global t
    
  print('time: '+ "{:.2f}".format(t/1e-9)+' ns of '+"{:.0f}".format(simulationtime/1e-9)+' ns ')
  time = 0
  while time<dtplot:
     time += dt        
     t += dt

     # inner nodes 
     for i in range(2,Xsize):      
        if i < int(boundary1/dx):
            alpha1 = alpha*1/(RefractiveIndexAmbient)
            alpha2 = alpha*1/(RefractiveIndexAmbient)
            gamma0 = 0
        elif i == int(boundary1/dx):
            alpha1 = alpha*1/(RefractiveIndexAmbient)
            alpha2 = alpha*1/(RefractiveIndexPlug)
            gamma0 = 0 
        elif (i > int(boundary1/dx)) and (i < int(boundary2/dx)):
            alpha1 = alpha*1/(RefractiveIndexPlug)
            alpha2 = alpha*1/(RefractiveIndexPlug)
            gamma0 = 0
        elif i == int(boundary2/dx):
            alpha1 = alpha*1/(RefractiveIndexPlug);
            alpha2 = alpha*1/(RefractiveIndexMedium);
            gamma0 = 0
        elif i > int(boundary2/dx):
            alpha1 = alpha*1/(RefractiveIndexMedium);
            alpha2 = alpha*1/(RefractiveIndexMedium);
            gamma0 = gamma

        line[i,2] = 1/(1+0.5*gamma0*dt)*(
                          2*line[i,1] -line[i,0] +
                          gamma0*dt/2*line[i,0]
                                     + alpha2**2*(line[i+1,1]-line[i,1])
                                     - alpha1**2*(line[i,1]-line[i-1,1])
                                       )


    # outer nodes 

     if BoundaryConditionLeft=='fixed':
          line[1,2] = 0
     elif BoundaryConditionLeft=='open':
          alpha1 = alpha*1/(RefractiveIndexAmbient)
          line[1,2] = -line[1,0] + 2*line[1,1] + alpha1**2*(line[2,1]-line[1,1]);
     elif BoundaryConditionLeft=='absorbing':
          alpha1 = alpha*1/(RefractiveIndexAmbient);
          line[1,2] = line[2,1] + (alpha1-1)/(alpha1+1)*(line[2,2]-line[1,1]);
 
     if BoundaryConditionRight=='fixed':
          line[Xsize,2] = 0
     elif BoundaryConditionRight=='open':
          alpha2 = alpha*1/(RefractiveIndexMedium);
          line[Xsize,2] = -line[Xsize,0] + 2*line[Xsize,1] - alpha2**2*(line[Xsize,1]-line[Xsize-1,1]);
     elif BoundaryConditionRight=='absorbing':
          alpha2 = alpha*1/(RefractiveIndexMedium);
          line[Xsize,2] = line[Xsize-1,1] + (alpha2-1)/(alpha2+1)*(line[Xsize-1,2]-line[Xsize,1])

     # write array back 
     for i in range(Xsize+1):
        line[i,0] = line[i,1];
        line[i,1] = line[i,2];
 
       
      
     lineamp.set_ydata(line[:,0]) 
     fig.suptitle('Time: ' + "{:.0f}".format(t/1e-9)+' ns')
  
  time = 0 
  

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationfile,fps=25,dpi=300)

print('End Simulation')

