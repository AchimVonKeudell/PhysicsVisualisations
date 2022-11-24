#
# Etching of a solid
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ----------------------------------------------------------------------------
# Monte Carlo Simulation etching 
# ----------------------------------------------------------------------------
#model = 'directional'
model = 'isotrop'

if model == 'directional':
    cosexponent = 20
    animationfile = 'etchdirectional.gif'
if model == 'isotrop':
    cosexponent = 1
    animationfile = 'etchisotrop.gif'

limitcf = 2
chemetch = 0.5
dt = 0.1

NBParticles = 5e4
NBParticlesReplot = 1e2
NBBinsX = 200
NBBinsY = 200
ms = 1 # markersize

# ---------------------------------------------------------------
# Setup structure
# ---------------------------------------------------------------
x = np.zeros(NBBinsX)
y = np.zeros(NBBinsX)
bins = np.zeros([NBBinsX,NBBinsY])

# material
ytop = int(0.75*NBBinsY)
for i in range(NBBinsX):
    for j in range(ytop):
        bins[i,j] = 1

# mask
ylayer = ytop + 5
for i in range(NBBinsX):
    for j in range(ytop,ylayer):
        bins[i,j] = 2

# large trench
xwidth = 50
xgaplarge = 0.7*NBBinsX
xgap = int(xgaplarge-xwidth/2)
for i in range(xgap,xgap+xwidth):
    for j in range(ytop,ylayer):
        bins[i,j] = 1

# small trench
xwidth = 20
xgapsmall = 0.3*NBBinsX
xgap = int(xgapsmall-xwidth/2)
for i in range(xgap,xgap+xwidth):
    for j in range(ytop,ylayer):
        bins[i,j] = 1

# ---------------------------------------------------------------
# Count surrounding CF species
# ---------------------------------------------------------------                
def surrounding(i,j):
    pt = 0
    if i>1 and i<NBBinsX-1 and j>1 and j<NBBinsY-1:    
      if bins[i-1,j-1] == 3: pt +=1
      if bins[i,j-1] == 3: pt +=1
      if bins[i+1,j-1] == 3: pt +=1
      if bins[i-1,j] == 3: pt +=1
      if bins[i+1,j] == 3: pt +=1
      if bins[i-1,j+1] == 3: pt +=1
      if bins[i,j+1] == 3: pt +=1
      if bins[i-1,j+1] == 3: pt +=1
    return pt

# ---------------------------------------------------------------
# center of mass surroudning species
# ---------------------------------------------------------------
def cmsurrounding(i,j):
    pt = 0
    xcm = 0
    ycm = 0
    if i>1 and i<NBBinsX-1 and j>1 and j<NBBinsY-1:    
      if bins[i-1,j-1] != 0: 
          xcm += 1*(-1)
          ycm += 1*(-1)
      if bins[i,j-1] != 0:
          xcm += 1*(0)
          ycm += 1*(-1)
      if bins[i+1,j-1] != 0:
          xcm += 1*(1)
          ycm += 1*(-1)
      if bins[i-1,j] != 0: 
          xcm += 1*(-1)
          ycm += 1*(0)
      if bins[i+1,j] != 0: 
          xcm += 1*(1)
          ycm += 1*(0)
      if bins[i-1,j+1] != 0: 
          xcm += 1*(-1)
          ycm += 1*(0)
      if bins[i+1,j+1] != 0:
          xcm += 1*(1)
          ycm += 1*(1)
          pt +=1
      if bins[i-1,j+1] != 0: 
          xcm += 1*(-1)
          ycm += 1*(1)
    return [xcm/8,ycm/8]

# ---------------------------------------------------------------
# center of mass surroudning species
# ---------------------------------------------------------------
def incidenceangle(x,y,xcm,ycm,vx,vy):
    return 0
        
# ---------------------------------------------------------------
# Setup Plot
# ---------------------------------------------------------------
# structure
fig, axp = plt.subplots(1,1,figsize=(4.5,4.5))
line1 = axp.plot(x,y,color='darkgray',marker='s',label='material',markersize=ms,lw=0)[0]
line2 = axp.plot(x,y,color='teal',marker='s',label='mask',markersize=ms,lw=0)[0]
line3 = axp.plot(x,y,color='orange',marker='s',label='CF',markersize=ms,lw=0)[0]
axp.set(xlabel="x",ylabel="z")
axp.set(xlim=(0,NBBinsX),ylim=(0,NBBinsY))

# initial distributions
#angle = np.linspace(0,np.pi,100)
#length = 20*np.cos(-np.pi/2+angle)**cosexponent
#xd = NBBinsX*0.7 + length*np.cos(angle)
#xd2 = NBBinsX*0.3 + length*np.cos(angle)
#yd = NBBinsY - length*np.sin(angle)
#linec1 = axp.plot(xd,yd,color='r',label='cos^n dist.',lw=1)        
#linec2 = axp.plot(xd2,yd,color='r',lw=1)        


# initial distributions ions
angle = np.linspace(0,np.pi,24)
length = 0.01+40*np.cos(-np.pi/2+angle)**cosexponent
xd = np.cos(angle) 
yd = -np.sin(angle) 
# one arrow for legend
axp.quiver(NBBinsX*0.7, NBBinsY, 0, -1, pivot = 'tail', scale_units='y',scale=1,color="r",label='ions', width=3e-3, alpha=0.5)
for i in range(len(xd)):
     axp.quiver(NBBinsX*0.7, NBBinsY, xd[i], yd[i], pivot = 'tail', scale_units='y',scale=1/length[i],color="r", width=3e-3, alpha=0.5)
     axp.quiver(NBBinsX*0.3, NBBinsY, xd[i], yd[i], pivot = 'tail', scale_units='y',scale=1/length[i],color="r", width=3e-3, alpha=0.5)

# initial distributions neutrals
angle = np.linspace(0,np.pi,12)
length = 0.01+20*np.cos(-np.pi/2+angle)
xd = np.cos(angle) #NBBinsX*0.5 +
yd = -np.sin(angle) # NBBinsY
axp.quiver(NBBinsX*0.7, NBBinsY, 0, -1, pivot = 'tail', scale_units='y',scale=1,color="orange",label='neutrals', width=3e-3,alpha=0.5)
for i in range(len(xd)):
     axp.quiver(NBBinsX*0.7, NBBinsY, xd[i], yd[i], pivot = 'tail', scale_units='y',scale=1/length[i],color="orange", width=3e-3,alpha=0.5)
     axp.quiver(NBBinsX*0.3, NBBinsY, xd[i], yd[i], pivot = 'tail', scale_units='y',scale=1/length[i],color="orange", width=3e-3,alpha=0.5)

axp.legend(loc=1,fontsize=6)


# setup info box
infobox = ''
infobox += 'cos_exponent: ' + "{:.0f}".format(cosexponent) + ' ' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
axp.text(0.02,0.98,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=axp.transAxes)

# ---------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------
# Initialize Simulation
Particles = 0
print('Start Simulation')

def animate(k):
    
    global bins,NBParticlesReplot
    
    # --------------------------------------------------------------------
    # Loop ions
    # --------------------------------------------------------------------
    NBP = 0
    while NBP<NBParticlesReplot:
        NBP += 1
        y0 = NBBinsY
        x0 = np.random.rand()*NBBinsX
        xp = x0
        yp = y0
        v0 = 1
        stheta = np.sqrt(np.random.rand())**cosexponent
        ctheta = np.sqrt(1-stheta*stheta)
        phi = 2*np.pi*np.random.rand()
        vy = -v0*ctheta
        vx = v0*stheta*np.cos(phi)
    
        t = 0 # start at t=0
        a = 0 # start without Force
        while yp>0:
            t += dt
            xp += vx*dt + 1/2*dt**2*a # new position        
            yp += vy*dt + 1/2*dt**2*a # new position 
    
            if xp<0:
               xp = NBBinsX+xp
            if xp>NBBinsX:
               xp = xp-NBBinsX
            xi = int(np.trunc(xp)) # x index bin
            yi = int(np.trunc(yp))  # y index bin
            if xi>=0 and xi<NBBinsX and yi>=0 and yi<NBBinsY:
               if bins[xi,yi] == 1 or bins[xi,yi] == 3: # hits material or etch layer
                   bins[int(np.trunc(xp)),int(np.trunc(yp))] = 0
                   break
               elif bins[xi,yi] == 2: # hits mask                  
                   break
               
    # --------------------------------------------------------------------
    # Loop radicals
    # --------------------------------------------------------------------
    NBP = 0
    while NBP<NBParticlesReplot:
        NBP += 1
        y0 = NBBinsY
        x0 = np.random.rand()*NBBinsX
        xp = x0
        yp = y0
        v0 = 1
        stheta = np.sqrt(np.random.rand())
        ctheta = np.sqrt(1-stheta*stheta)
        phi = 2*np.pi*np.random.rand()
        vy = -v0*ctheta
        vx = v0*stheta*np.cos(phi)
    
        t = 0 # start at t=0
        a = 0 # start without Force
        while yp>0:
            t += dt
            xp += vx*dt + 1/2*dt**2*a # new position        
            yp += vy*dt + 1/2*dt**2*a # new position 
    
            if xp<0:
               xp = NBBinsX+xp
            if xp>NBBinsX:
               xp = xp-NBBinsX
            xi = int(np.trunc(xp)) # x index bin
            yi = int(np.trunc(yp))  # y index bin
            if xi>=1 and xi<NBBinsX and yi>=0 and yi<NBBinsY:
               if bins[xi,yi] != 0: 
                   if surrounding(xi,yi)<limitcf:
                      # move a step back 
                      xi = int(xp - vx*dt*0.5)
                      yi = int(yp - vy*dt*0.5)
                      if xi<0: xi = NBBinsX+xi
                      if xi>NBBinsX-1: xi = xi-(NBBinsX-1)           
                      if yi<0: yi = NBBinsY+yi
                      if yi>NBBinsX-1: yi = yi-(NBBinsY-1) 
                      # add radical only if place empty and maximal 2 on mask
                      if bins[xi,yi] == 0 and yi<ytop+2: bins[xi,yi] = 3
                   break
                   
    print('Step: ',int(k*NBParticlesReplot),' of ', int(NBParticles))        
    # material
    x = []
    y = []
    for i in range(NBBinsX):
        for j in range(NBBinsY):
            if bins[i,j] == 1:
               x.append(i+0.5)
               y.append(j+0.5)               
    line1.set_xdata(x)
    line1.set_ydata(y) 
    
    # mask
    xl = []
    yl = []
    for i in range(NBBinsX):
        for j in range(NBBinsY):
            if bins[i,j] == 2:
               xl.append(i+0.5)
               yl.append(j+0.5)               
    line2.set_xdata(xl)
    line2.set_ydata(yl)     
 
    # radicals
    x2 = []
    y2 = []
    for i in range(NBBinsX):
        for j in range(NBBinsY):
            if bins[i,j] == 3:
               x2.append(i+0.5)
               y2.append(j+0.5)               
    line3.set_xdata(x2)
    line3.set_ydata(y2)     
     
    fig.suptitle('Particles: ' + "{:.0f}".format(int(k*NBParticlesReplot)))
    
anim = animation.FuncAnimation(fig,animate,interval=1,frames=int(NBParticles/NBParticlesReplot))
anim.save(animationfile,fps=25,dpi=180)


x = []
y = []
for i in range(NBBinsX):
    for j in range(NBBinsY):
            if bins[i,j] == 1:
               x.append(i+0.5)
               y.append(j+0.5)  
xl = []
yl = []
for i in range(NBBinsX):
        for j in range(NBBinsY):
            if bins[i,j] == 2:
               xl.append(i+0.5)
               yl.append(j+0.5)    
               
x2 = []
y2 = []
for i in range(NBBinsX):
        for j in range(NBBinsY):
            if bins[i,j] == 3:
               x2.append(i+0.5)
               y2.append(j+0.5)               
                
#axp.plot(x,y,color='darkgray',marker='s',markersize=ms,lw=0)[0]    
#axp.plot(xl,yl,color='teal',marker='s',markersize=ms,lw=0)[0]    
#axp.plot(x2,y2,color='orange',marker='s',markersize=ms,lw=0)[0]    
