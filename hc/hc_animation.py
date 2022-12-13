#
# Hollow cathode effect
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

model = 'hc DC'

xmax = 100
ymax = 100
dx = 1
dy = 1

steps = 10000
animationname = 'hc.mp4'
dt = 0.1
dtreplot = 1

tmax = 500

v0 = 0.001
qdm = 10000
width = 0.4*xmax
depth = 0.7*ymax

NBparticles = 10

potential = np.zeros((xmax,ymax))
potentialimage = np.zeros((xmax,ymax))
efeldx = np.zeros((xmax,ymax))
efeldy = np.zeros((xmax,ymax))
efeldximage = np.zeros((xmax,ymax))
efeldyimage = np.zeros((xmax,ymax))
potentialnew = np.zeros((xmax,ymax))

# --------------------------------------
#  Solve Poisson Equation
#   Gauss Seidel method
# --------------------------------------
def SolvePoisson(k):
    global efeldx,efeldy,potential
    poissoniteration = 0
    poissoniterationlimit = 1000
    while poissoniteration<poissoniterationlimit:
        poissoniteration += 1
        # Gauss Seidel step
        for i in range(2,xmax-1):
            for j in range(2,ymax-1):
                potentialnew[i,j] = 0.25*(potential[i-1,j]+potential[i+1,j]
                                          +potential[i,j-1]+potential[i,j+1])
        
        # set boundary
        for i in range(xmax):
           for j in range(ymax):
              if ((i< int(xmax/2-width/2)) or (i> int(xmax/2+width/2))):
                 if j<int(depth):
                    potentialnew[i,j] = 1
              if j == 0 or j == 1: potentialnew[i,j] = 1         

        # set central part
        for i in range(xmax):
           for j in range(ymax):
              if ((i> int(xmax/2-width/7)) and (i< int(xmax/2+width/7))) and j>0.1*ymax:
                    potentialnew[i,j] = 0
              if j == 0 or j == 1: potentialnew[i,j] = 1         

                
        # write array back and check for deviation
        diffmax = 0;
        for i in range(xmax):
            for j in range(ymax):
               if abs(potential[i,j]-potentialnew[i,j])>diffmax:
                   diffmax = abs(potential[i,j]-potentialnew[i,j])
               potential[i,j] = potentialnew[i,j]
        if diffmax < 1/xmax**2*10: return # return when diff < threshold
        if k == poissoniterationlimit: return # return when max iter. reached


# --------------------------------------
#  Generate Structure, fields
# --------------------------------------
for i in range(xmax):
    for j in range(ymax):
        if ((i< int(xmax/2-width/2)) or (i> int(xmax/2+width/2))):
            if j<int(depth):
                potential[i,j] = 1
        if j == 0: potential[i,j] = 1         

print('Solve Poisson Equation')
SolvePoisson(1) 
potential = 1 - potential
for i in range(2,xmax-1):
  for j in range(2,ymax-1):
      efeldx[i,j] = -(potential[i+1,j]-potential[i-1,j])/2
      efeldy[i,j] = -(potential[i,j+1]-potential[i,j-1])/2

potentialimage = np.transpose(potential)

# --------------------------------------
#  PIC, fields to force
# --------------------------------------
def force(r):
    global efeldx,efeldy,xmax
   # if r[0]>xmax/2-width/2 and r[0]<xmax/2+width/2:
    Ex = (r[0]-xmax/2)**3*1e-7
   # else:
   #     Ex = 0
    #print(Ex)
    return np.array([Ex,0])

def forcepic(r):
    global efeldx,efeldy
        
    ibin = int(r[0]/dx)
    ibinp1 = ibin+1
    jbin = int(r[1]/dy)
    jbinp1 = jbin+1
    
    if ibin>=0 and ibin<xmax-1 and jbin>=0 and jbin<ymax-1:
      xi = ibin*dx
      xip1 = ibin*dx + dx
      yj = jbin*dy 
      yjp1 = jbin*dy + dy

      A4 = (xip1-r[0])*(yjp1-r[1])/(dx*dy)
      A3 = (r[0]-xi)*(yjp1-r[1])/(dx*dy)
      A2 = (xip1-r[0])*(r[1]-yj)/(dx*dy)
      A1 = (r[0]-xi)*(r[1]-yj)/(dx*dy)

      Ex = (A4*efeldx[ibin,jbin]+
          A3*efeldx[ibinp1,jbin]+
          A2*efeldx[ibin,jbinp1]+
          A1*efeldx[ibinp1,jbinp1])
      Ey = (A4*efeldy[ibin,jbin]+
          A3*efeldy[ibinp1,jbin]+
          A2*efeldy[ibin,jbinp1]+
          A1*efeldy[ibinp1,jbinp1])    
    else:
      Ex = 0
      Ey = 0

    return np.array([Ex,Ey])

 

#----------------------------------
# Define Plot
#----------------------------------
fig = plt.figure(constrained_layout=True)
ax = []
ax.append(fig.add_subplot())

# potential
ax[0].set(xlim=(0,xmax),ylim=(0,ymax),xlabel='x',ylabel='y')
im = plt.imshow(potentialimage,cmap='plasma',label='potential')
fig.colorbar(im,ax=ax[0],label='potential')

# trajectories
linea = []
lineap = []
xt =[]
yt = []
for i in range(NBparticles):
  linea.append(ax[0].plot(xt,yt,color='r',lw=1, label ='')[0])
  lineap.append(ax[0].plot(xt,yt,color='r',lw=0, label ='',marker='o',markersize=2)[0])

# contour trench
xt = [0,xmax/2-width/2,xmax/2-width/2,xmax/2+width/2,xmax/2+width/2,xmax]
yt = [depth,depth,1,1,depth,depth] 
ax[0].plot(xt,yt,color='w',lw=2)
#ax[0].legend(fontsize=6)

# info box
infobox = ''
infobox += 'model: ' + model + '' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8) 
ax[0].text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax[0].transAxes)


# ----------------------------------
# Main loop
# ----------------------------------
t = 0
time = 0

trajectorypt = 0      

angle = np.random.rand()*np.pi/2*0.3+np.pi/2*0.7
vx = v0*np.cos(angle)
if np.random.rand()>0.5:
   vy = v0*np.sin(angle)
else:   
   vy = -v0*np.sin(angle)

if np.random.rand()>0.5:
    x = xmax/2-width/2+0.1
    y = np.random.rand()*depth+1
    r = np.array([x,y])
else:   
    x = xmax/2+width/2-0.1
    y = np.random.rand()*depth+1
    r = np.array([x,y])
    vx = -vx


v = np.array([vx,vy])   
a = -forcepic(r)
dr = np.array([0,0])
    
xt = []
yt = []

print('Start Particle 0')
   

def animate(k):

    global t,time,dt,r,v,a,xt,yt,dr
    global trajectorypt
    # Gauss Seidel solution
    
    if trajectorypt>NBparticles-1: return
    
    if (r[0]<xmax and r[0]>0 and r[1]<ymax and r[1]>0) and t<tmax:
      if t == 0: 
          print('r begin: ',r)  
          print('v begin: ',v)  
          print('a begin: ',a)  
      while time<dtreplot:
        t += dt
        time += dt  
      
        dr = dt*v+0.5*dt**2*a
        r = r + dr
        v = v + 0.5*dt*a
        a = -forcepic(r)
        v = v + 0.5*dt*a
      time = 0  
      if r[0]>0 or r[0]<xmax or r[1]>0 or r[1]<ymax: 
        xt.append(r[0])
        yt.append(r[1])

      linea[trajectorypt].set_xdata(xt)
      linea[trajectorypt].set_ydata(yt)
      lineap[trajectorypt].set_xdata(r[0])
      lineap[trajectorypt].set_ydata(r[1])      
    else: # new particles
      linea[trajectorypt].set_color('b')
      linea[trajectorypt].set_linestyle('dashed')
      lineap[trajectorypt].set_color('b')
      
      trajectorypt += 1  
      print('Start Particle ',trajectorypt)
      #print('time last particle ',t)
      
      t = 0      
      time = 0   
      angle = np.random.rand()*np.pi/2*0.3+np.pi/2*0.7
      vx = v0*np.cos(angle)
      if np.random.rand()>0.5:
         vy = v0*np.sin(angle)
      else:   
         vy = -v0*np.sin(angle)

      if np.random.rand()>0.5:
         x = xmax/2-width/2+0.1
         y = np.random.rand()*depth+1
         r = np.array([x,y])
      else:   
         x = xmax/2+width/2-0.1
         y = np.random.rand()*depth+1
         r = np.array([x,y])
         vx = -vx
  
      v = np.array([vx,vy])
    
      a = -forcepic(r)
      dr = np.array([0,0])
    
      #print('r0:',r)
      #print('v0:',v)
      xt = []
      yt = []
   
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)