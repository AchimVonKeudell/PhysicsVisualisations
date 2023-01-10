#
# Waves 2d
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
model = 'open boundaries'
model = 'traveling peak'
model = 'plane wave'
model = 'diffractionslit'
model = 'diffractiondoubleslit'

dtplot = 2e-4 # time distance in between frames for animation
simulationtime = 400e-3 # total length of simulation in s
dt=1e-4 # time constant for iteration
gamma = 0 # damping in the medium
omega = 1e3 # frequency of plane wave
amplitude = 0.05
    
if model == 'fixed boundaries':
    animationfile = 'wave2d_fixed.mp4'
    BoundaryCondition = 'fixed'
    InitialCondition = 'Zero'
if model == 'open boundaries':
    animationfile = 'wave2d_open.mp4'
    BoundaryCondition = 'open'
    InitialCondition = 'Zero'
if model == 'absorbing boundaries':
    animationfile = 'wave2d_absorbing.mp4'
    BoundaryCondition = 'absorbing'
    InitialCondition = 'Zero'
if model == 'traveling peak':
    animationfile = 'wave2d_travelingpeak.mp4'
    BoundaryCondition = 'absorbing'
    InitialCondition = 'Direction'
if model == 'plane wave':
    animationfile = 'wave2d_planewave.mp4'
    BoundaryCondition = 'absorbing'
    InitialCondition = 'Line'
    Slit = False
if model == 'diffractionslit':
    animationfile = 'wave2d_slit.mp4'
    BoundaryCondition = 'absorbing'
    InitialCondition = 'Line'
    Slit = True
    DoubleSlit = False
if model == 'diffractiondoubleslit':
    animationfile = 'wave2d_doubleslit.mp4'
    BoundaryCondition = 'absorbing'
    InitialCondition = 'Line'
    Slit = False
    DoubleSlit = True
    
# --------------------------------------------------------------------------
# Define initial constants
# --------------------------------------------------------------------------
time = 0
t = 0
steps = int(simulationtime/dtplot) # total number of frames

Xsize = 100
peakwidth = 3 # in pixels
dx = 0.2
SlitPositionX = 20
SlitPositionY = 50
SlitWidth = 5
DoubleSlitDistance = 30

field = np.zeros([Xsize+1,Xsize+1,3])
dfield = np.zeros([Xsize+1,Xsize+1])
x = np.linspace(0,dx*Xsize,Xsize+1) 
X, Y = np.meshgrid(x, x)

# --------------------------------------------------------------------------
# Initialize Distributions
# The initial distribution needs to be electrical neutral to
# allow the calculation of the electrical potential for periodic
# boundary condition
# --------------------------------------------------------------------------

# Starting Condition Amplitude 
if InitialCondition == 'Zero':
  for i in range(Xsize+1): 
    for j in range(Xsize+1): 
      field[i,j,0] = np.exp(-((i-Xsize/2)**2+(j-Xsize/2)**2)/(peakwidth/dx)**2)
if InitialCondition == 'Direction':
  for i in range(Xsize+1): 
    for j in range(Xsize+1): 
      field[i,j,0] = np.exp(-((i-Xsize/2)**2+(j-Xsize/2)**2)/(peakwidth/dx)**2)
if InitialCondition == 'Line':
  for i in range(Xsize+1): 
    for j in range(Xsize+1): 
      field[i,j,0] = 0

   
vphase = 3e2
gamma0 = 0
alpha = vphase*dt/dx
if alpha>1: print('von Neumann condition violated')
 
    
# Starting Condition Velocity 
# runs with vphase to x,y
for i in range(2,Xsize-1):
  for j in range(2,Xsize-1):
      if InitialCondition == 'Zero':
         dfield[i,j] = 0
      elif InitialCondition == 'Direction':
         dfield[i,j] = -vphase*((field[i+1,j,0]-field[i-1,j,0])/(2*dx)) 
      elif InitialCondition == 'Line':
         dfield[i,1] = 0.02*np.cos(0*omega)*omega  
         #dfield[i,1] = -vphase*((field[i+1,1,0]-field[i-1,1,0])/(2*dx)) 

# Starting Condition Ghost Point in Time }
for i in range(Xsize-1):
  for j in range(Xsize-1):
      field[i,j,1] = field[i,j,0]+0.5*alpha**2*((field[i+1,j,0]-2*field[i,j,0]+field[i-1,j,0])+
                                                (field[i,j+1,0]-2*field[i,j,0]+field[i,j-1,0]))+ dfield[i,j]*dt
  
 
# --------------------------------------------------------------------------
# Define plots
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4.5))
ax = fig.add_subplot(111, projection='3d')

plot = [ax.plot_surface(X, Y, field[:,:,0], color='0.75', rstride=1, cstride=1,cmap="magma",vmin=-1,vmax=1)]
ax.set(xlabel='x (m)',ylabel='y (m)')
ax.set_zlim(-1,1)
ax.set_zticks([])

ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

#im = plt.imshow(field,extent=(0,Xsize*dx,0,Xsize*dx),cmap='RdBu',label='abs(E)',alpha=0.5)
fig.colorbar(plot[0],ax=ax,label='amplitude',cmap="magma",extend='both')

# info box
infobox = ''
infobox += 'model: ' + model +'\n'
infobox += 'v_phase: '+ "{:.0f}".format(vphase)+' (m/s)'
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text2D(0.02,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)


print('Start Simulation of '+"{:.0f}".format(steps)+' frames')

# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------
t = 0
def animate(k):
  global t
    
  print('time: '+ "{:.2f}".format(t/1e-3)+' ms of '+"{:.0f}".format(simulationtime/1e-3)+' ms ')
  time = 0
  while time<dtplot:
     time += dt        
     t += dt

     # inner nodes 
     for i in range(2,Xsize):      
      for j in range(2,Xsize):      
        gamma0 = gamma
        field[i,j,2] = 1/(1+0.5*gamma0*dt)*(
                          2*field[i,j,1] -field[i,j,0] + gamma0*dt/2*field[i,j,0]
                                     + alpha**2*(field[i+1,j,1]-field[i,j,1])
                                     - alpha**2*(field[i,j,1]-field[i-1,j,1])
                                     + alpha**2*(field[i,j+1,1]-field[i,j,1])
                                     - alpha**2*(field[i,j,1]-field[i,j-1,1])
                                       )

    # outer nodes 
     if BoundaryCondition=='fixed':
        for i in range(Xsize+1):
          field[i,1,2] = 0
          field[i,Xsize,2] = 0
          field[1,i,2] = 0
          field[Xsize,i,2] = 0
     elif BoundaryCondition=='open':
        for i in range(Xsize+1):
          field[1,i,2] = -field[1,i,0] + 2*field[1,i,1] + alpha**2*(field[2,i,1]-field[1,i,1]);
          field[Xsize,i,2] = -field[Xsize,i,0] + 2*field[Xsize,i,1] - alpha**2*(field[Xsize,i,1]-field[Xsize-1,i,1])
          field[i,1,2] = -field[i,1,0] + 2*field[i,1,1] + alpha**2*(field[i,2,1]-field[i,1,1]);
          field[i,Xsize,2] = -field[i,Xsize,0] + 2*field[i,Xsize,1] - alpha**2*(field[i,Xsize,1]-field[i,Xsize-1,1])
     elif BoundaryCondition=='absorbing':
        for i in range(Xsize+1):    
          field[1,i,2] = field[2,i,1] + (alpha-1)/(alpha+1)*(field[2,i,2]-field[1,i,1])
          field[i,1,2] = field[i,2,1] + (alpha-1)/(alpha+1)*(field[i,2,2]-field[i,1,1])
          field[Xsize,i,2] = field[Xsize-1,i,1] + (alpha-1)/(alpha+1)*(field[Xsize-1,i,2]-field[Xsize,i,1])
          field[i,Xsize,2] = field[i,Xsize-1,1] + (alpha-1)/(alpha+1)*(field[i,Xsize-1,2]-field[i,Xsize,1])

     # excitation 
     if InitialCondition == 'Line':
       for i in range(1,Xsize+1):      
          field[i,1,2] = amplitude*np.sin(t*omega)
   
     if Slit:
       for i in range(1,Xsize+1):
           if i<SlitPositionY-0.5*SlitWidth:
             field[i,SlitPositionX,2] = 0
           if i>SlitPositionY+0.5*SlitWidth:
             field[i,SlitPositionX,2] = 0
     
     if DoubleSlit:
       for i in range(1,Xsize+1):
           if i<SlitPositionY-DoubleSlitDistance*0.5-0.5*SlitWidth:
             field[i,SlitPositionX,2] = 0
           if (i>SlitPositionY-DoubleSlitDistance*0.5+0.5*SlitWidth
             and i<SlitPositionY+DoubleSlitDistance*0.5-0.5*SlitWidth):
             field[i,SlitPositionX,2] = 0
           if i>SlitPositionY+DoubleSlitDistance*0.5+0.5*SlitWidth:
             field[i,SlitPositionX,2] = 0
             
         

     # write array back 
     for i in range(Xsize+1):
      for j in range(Xsize+1):
        field[i,j,0] = field[i,j,1];
        field[i,j,1] = field[i,j,2];
 
     plot[0].remove()
     plot[0] = ax.plot_surface(X[1:,1:], Y[1:,1:], field[1:,1:,0], cmap="magma")
     #im.set_array(np.transpose(field[:,:,0]))  
      
#     lineamp.set_ydata(line[:,0]) 
     fig.suptitle('Time: ' + "{:.0f}".format(t/1e-3)+' ms')
  
  time = 0 
  

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationfile,fps=25,dpi=300)

print('End Simulation')

