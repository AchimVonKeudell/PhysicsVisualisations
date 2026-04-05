#
# 1d wave on a line
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
import matplotlib.pyplot as plt # for plot functons
import scipy as sp
import matplotlib.animation as animation


model = 'pulse'
model = 'pulsewidths'
#model = 'signal'
#model = 'signalvphase'
#model = 'signalhighfrequency'

dtplot = 2e-10 # time distance in between frames for animation
simulationtime = 400e-9 # total length of simulation in s
dt=1e-12 # time constant for iteration
omega = 3e8*np.pi
f = omega/(2*np.pi)
period = 1/f
pulsewidth1 = 0.1
pulsewidth2 = 0.1

RefractiveIndexMedium = 1
gamma = 0 # damping in the medium
    
if model == 'pulse':
    animationfile = 'line1d_pulse.mp4'
    ffactor = 1
    zoomfactor = 1
    vphasefactor = 1
    boundary1 = 'open'
    boundary2 = 'fixed'
    showarrow = False
if model == 'pulsewidths':
    animationfile = 'line1d_pulsewidth.mp4'
    ffactor = 1
    zoomfactor = 1
    vphasefactor = 1
    boundary1 = 'open'
    boundary2 = 'open'
    showarrow = True
    omega = 3e8
    f = omega/(2*np.pi)
    period = 1/f
    pulsewidth1 = 0.1
    pulsewidth2 = 0.19
if model == 'signal':
    animationfile = 'line1d_signal.mp4'
    ffactor = 1
    zoomfactor = 1
    vphasefactor = 1 
    boundary1 = 'open'
    boundary2 = 'fixed'
    showarrow = False
if model == 'signalvphase':
    animationfile = 'line1d_signalvphase.mp4'
    ffactor = 1
    zoomfactor = 1
    vphasefactor = 2 
    boundary1 = 'open'
    boundary2 = 'open'
    showarrow = False    
if model == 'signalhighfrequency':
    animationfile = 'line1d_signalhigh.mp4'
    ffactor = 5
    zoomfactor = 1
    vphasefactor = 1 
    boundary1 = 'open'
    boundary2 = 'fixed'
    showarrow = False
    
# --------------------------------------------------------------------------
# Define initial constants
# --------------------------------------------------------------------------
time = 0
t = 0
steps = int(simulationtime/dtplot) # total number of frames

#Xsize = 10000
#dx = 0.002

Xsize = 1000
dx = 0.02
x = np.linspace(0,dx*Xsize,Xsize+1)

# vector open
line = np.zeros([Xsize+1,3])

# vector fixed
line1 = np.zeros([Xsize+1,3])
  

# --------------------------------------------------------------------------
# Initialize Distributions
# The initial distribution needs to be electrical neutral to
# allow the calculation of the electrical potential for periodic
# boundary condition
# --------------------------------------------------------------------------

# Starting Condition Amplitude 
for i in range(Xsize): 
    line[i,0] = 0
    line1[i,0] = 0
   
vphase = 3e8
gamma0 = 0
alpha = vphase*dt/dx

if alpha>1: print('von Neumann condition violated')
  
  
# Starting Condition Ghost Point in Time }
for i in range(Xsize-1):
      line[i,1] = line[i,0]+0.5*alpha**2*(line[i+1,0]-2*line[i,0]+line[i-1,0]);
#      line1[i,1] = line1[i,0]+0.5*alpha**2*(line1[i+1,0]-2*line1[i,0]+line1[i-1,0]);
  
 
# --------------------------------------------------------------------------
# Define plots
# --------------------------------------------------------------------------
fig, ax = plt.subplots(1,1,figsize=(8,4.5))
plt.subplots_adjust(right=0.85)
yrange = 4*1/zoomfactor
yoffset = 1.5*1/zoomfactor

ax.set(xlabel='position (m)',ylabel='amplitude',ylim=(-yrange,yrange),xlim=(-0.1*Xsize*dx,1.1*Xsize*dx))

# line open boundary
lineopen = ax.plot(x,line[:,0],color='b')[0]
lineopena = ax.plot([0,0],[yoffset,yoffset],color='r',linestyle='dashed')[0]
lineopenap = ax.plot(0,0,color='r',marker='o')[0]
if showarrow:
    lineopenquiver = ax.quiver(x[int(Xsize/2)],yoffset,0,1,color='g',label='surfer')
if boundary1 == 'open':
   lineopene = ax.plot(Xsize*dx,0,color='b',marker='o',fillstyle='none')[0]
else:
   lineopene = ax.plot(Xsize*dx,0,color='b',marker='o',fillstyle='full')[0]

# line fixed boundary
linefixed = ax.plot(x,line1[:,0],color='b')[0]
linefixeda = ax.plot([0,0],[-yoffset,-yoffset],color='r',linestyle='dashed')[0]
linefixedap = ax.plot(0,0,color='r',marker='o')[0]
if showarrow:
    linefixedquiver = ax.quiver(x[int(Xsize/2)],-yoffset,0,1,color='g')
if boundary2 == 'open':
    linefixede = ax.plot(Xsize*dx,0,color='b',marker='o',fillstyle='none')[0]
else:    
    linefixede = ax.plot(Xsize*dx,0,color='b',marker='o',fillstyle='full')[0]


ax.legend(loc=1,fontsize=8)

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
infobox += 'model: ' + model +'\n'
infobox += 'boundary: ' + boundary1 + '\n'
infobox += 'vphase: ' + "{:.1f}".format(1) + ' x c \n'
infobox += 'f: ' + "{:.2e}".format(f)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)

infobox = ''
infobox += 'model: ' + model +'\n'
infobox += 'boundary: ' + boundary2 + '\n'
infobox += 'vphase: ' + "{:.1f}".format(1/vphasefactor) + ' x c \n'
infobox += 'f: ' + "{:.2e}".format(f)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.5,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)


print('Start Simulation of '+"{:.0f}".format(steps)+' frames')

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
     
     # excitation left side
     if model == 'signal' or model == 'signalhighfrequency' or model == 'signalvphase':
       amp = np.sin(ffactor*omega*t)
       line[1,0] = amp
       line[1,1] = amp
       line[1,2] = amp
       line1[1,0] = amp
       line1[1,1] = amp
       line1[1,2] = amp
     elif model == 'pulse' or model =='pulsewidths':    
       if t<=period:
         amp = np.exp(-(t-0.5*period)**2/(2*(period*pulsewidth1)**2))  
         line[1,0] = amp
         line[1,1] = amp
         line[1,2] = amp
         amp = np.exp(-(t-0.5*period)**2/(2*(period*pulsewidth2)**2))  
         line1[1,0] = amp
         line1[1,1] = amp
         line1[1,2] = amp

     # inner nodes 
     for i in range(2,Xsize):      
        alpha1 = alpha*1/(RefractiveIndexMedium);
        
        alpha2 = alpha*1/(RefractiveIndexMedium)*1/vphasefactor;
        gamma0 = gamma

        line[i,2] = 1/(1+0.5*gamma0*dt)*(
                          2*line[i,1] -line[i,0] +
                          gamma0*dt/2*line[i,0]
                                     + alpha1**2*(line[i+1,1]-line[i,1])
                                     - alpha1**2*(line[i,1]-line[i-1,1])
                                       )
        line1[i,2] = 1/(1+0.5*gamma0*dt)*(
                          2*line1[i,1] -line1[i,0] +
                          gamma0*dt/2*line1[i,0]
                                     + alpha2**2*(line1[i+1,1]-line1[i,1])
                                     - alpha2**2*(line1[i,1]-line1[i-1,1])
                                       )


     # outer nodes 
     #if BoundaryConditionLeft=='fixed':
     if boundary2 == 'fixed':    
       line1[Xsize,2] = 0
       line1[1,2] = 0
     else:  
       line1[1,2] = -line1[1,0] + 2*line1[1,1] + alpha2**2*(line1[2,1]-line1[1,1]);
       line1[Xsize,2] = -line1[Xsize,0] + 2*line1[Xsize,1] - alpha2**2*(line1[Xsize,1]-line1[Xsize-1,1]);
     #elif BoundaryConditionLeft=='open':
     if boundary1 == 'fixed':
       line[Xsize,2] = 0
       line[1,2] = 0
     else:  
       line[1,2] = -line[1,0] + 2*line[1,1] + alpha1**2*(line[2,1]-line[1,1]);
       line[Xsize,2] = -line[Xsize,0] + 2*line[Xsize,1] - alpha1**2*(line[Xsize,1]-line[Xsize-1,1]);
          
     #elif BoundaryConditionRight=='absorbing':
     #     line[Xsize,2] = line[Xsize-1,1] + (alpha2-1)/(alpha2+1)*(line[Xsize-1,2]-line[Xsize,1])

     # write array back 
     for i in range(Xsize+1):
        line[i,0] = line[i,1];
        line[i,1] = line[i,2];
        line1[i,0] = line1[i,1];
        line1[i,1] = line1[i,2];
        
     lineopen.set_ydata(yoffset+line[:,0]) 
     lineopena.set_ydata([yoffset,yoffset+line[1,0]]) 
     lineopenap.set_ydata(yoffset+line[1,0]) 
     lineopene.set_ydata(yoffset+line[Xsize,0])
     if showarrow:
       lineopenquiver.set_offsets([x[int(Xsize/2)],yoffset+line[int(Xsize/2),0]])
     
     linefixed.set_ydata(-yoffset+line1[:,0]) 
     linefixeda.set_ydata([-yoffset,-yoffset+line1[1,0]]) 
     linefixedap.set_ydata(-yoffset+line1[1,0]) 
     linefixede.set_ydata(-yoffset+line1[Xsize,0])
     if showarrow:
       linefixedquiver.set_offsets([x[int(Xsize/2)],-yoffset+line1[int(Xsize/2),0]])

     fig.suptitle('Time: ' + "{:.0f}".format(t/1e-9)+' ns')
  
  time = 0 
  

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationfile,fps=25,dpi=300)

print('End Simulation')

