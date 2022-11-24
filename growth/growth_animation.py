#
# Roughness evolution
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
#model = 'ballisticcosine'
#model = 'ballisticnormal'
model = 'surfacediffusion'
#model = 'isotrop'

if model == 'ballisticcosine':
    modeltext = 'Ballistic - cosine'
    cosexponent = 20
    theoryexponent = 0.5
    NBParticles = 1e4
    NBParticlesReplot = 1e2
    NBParticlesChangeColor =5e3
    animationfile = 'growth_ballisticcosine.gif'
if model == 'ballisticnormal':
    modeltext = 'Ballistic - normal'
    cosexponent = 0
    theoryexponent = 0.5
    NBParticles = 2e4
    NBParticlesReplot = 1e2
    NBParticlesChangeColor =5e3    
    animationfile = 'growth_ballisticnormal.gif'
if model == 'surfacediffusion':
    modeltext = 'Surface diffusion'
    cosexponent = 0
    theoryexponent = 0.25
    surfacesteps = 10
    NBParticles = 3e4
    NBParticlesReplot = 1e2
    NBParticlesChangeColor =4e3
    animationfile = 'growth_surfacediffusion_10.gif'


NBBinsX = 200
NBBinsY = 200
ms = 1 # markersize


# ---------------------------------------------------------------
# Setup structure
# ---------------------------------------------------------------
x = np.zeros(NBBinsX)
y = np.zeros(NBBinsX)
bins = np.zeros([NBBinsX,NBBinsY])

hprofilex = np.linspace(0,NBBinsX,NBBinsX)
hprofiley = np.ones(NBBinsX)

# material
ytop = 1
for i in range(NBBinsX):
    for j in range(ytop):
        bins[i,j] = 1

# ---------------------------------------------------------------
# center of mass surrounding species
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
# plot structure
fig, ax = plt.subplots(1,2,figsize=(8,4.5))
line1 = ax[0].plot(x,y,color='orange',marker='s',label='',markersize=ms,lw=0)[0]
line2 = ax[0].plot(x,y,color='olive',marker='s',label='',markersize=ms,lw=0)[0]
linehprofile = ax[0].plot(hprofilex,hprofiley,color='b',label='h',lw=1)[0]
ax[0].set(xlabel="x",ylabel="h")
ax[0].set(xlim=(0,NBBinsX),ylim=(0,NBBinsY))

# initial distributions
angle = np.linspace(0,np.pi,36)
if model == 'ballisticcosine':
  length = 0.01+40*np.cos(-np.pi/2+angle)**cosexponent
  xd = np.cos(angle) #NBBinsX*0.5 +
  yd = -np.sin(angle) # NBBinsY
  ax[0].quiver(NBBinsX/2, NBBinsY, 0, -1, pivot = 'tail', scale_units='y',scale=1/40,color="r", width=3e-3, label='incident species',alpha=0.5)
  for i in range(len(xd)):
     ax[0].quiver(NBBinsX/2, NBBinsY, xd[i], yd[i], pivot = 'tail', scale_units='y',scale=1/length[i],color="r", width=3e-3, alpha=0.5)
if model == 'ballisticnormal' or model == 'surfacediffusion':
  ax[0].quiver(NBBinsX/2, NBBinsY, 0, -1, pivot = 'tail', scale_units='y',scale=1/40,color="r", width=3e-3, label='incident species',alpha=0.5)
#linec1 = ax[0].plot(xd,yd,color='r',label='cos^n dist.',lw=1)        
ax[0].legend(loc=1,fontsize=6)

# setup info box
infobox = ''
infobox += modeltext + '\n' 
if model == 'surfacediffusion': infobox += 'surface steps: ' + "{:.0f}".format(surfacesteps) + '\n' 
infobox += 'cos_exponent: ' + "{:.0f}".format(cosexponent) + ' ' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[0].text(0.02,0.98,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax[0].transAxes)

# plot change height
ht = [0]
hmean = [0]
linehmean = ax[1].plot(ht,hmean,color='g',marker='s',label='sigma h',markersize=ms,lw=2)[0]
httheory = np.linspace(0,NBParticles,int(NBParticles))
hmeantheory = []
for i in range(int(NBParticles)):
    hmeantheory.append((i/NBBinsX)**theoryexponent)
linehmeantheroy = ax[1].plot(httheory,hmeantheory,color='black',label='theory',linestyle='dashed',lw=1)[0]
ax[1].set(xlabel="particles",ylabel="sigma h")
ax[1].set(xlim=(0,NBParticles),ylim=(0,NBBinsY**theoryexponent*2))
ax[1].legend(loc=1,fontsize=6)

# setup info box
infobox = ''
infobox += 'model: ' + modeltext + '\n' 
infobox += 'theory_exponent: ' + "{:.2f}".format(theoryexponent) + ' ' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[1].text(0.02,0.98,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax[1].transAxes)


# ---------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------
# Initialize Simulation
dt = 0.1
print('Start Simulation')

def animate(k):
    
    global bins,NBParticlesReplot,ColorPt,color
                  
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
        if model == 'ballisticcosine':
          stheta = np.sqrt(np.random.rand())
          ctheta = np.sqrt(1-stheta*stheta)
          phi = 2*np.pi*np.random.rand()
          vy = -v0*ctheta
          vx = v0*stheta*np.cos(phi)
        if model == 'ballisticnormal':
          vy = -v0
          vx = 0
        if model == 'surfacediffusion':
          vy = -v0
          vx = 0
    
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
               if bins[xi,yi] != 0: # hit
                      # move a step back 
                      xi = int(xp - vx*dt)
                      yi = int(yp - vy*dt)
                      if xi<0: xi = NBBinsX-1+xi
                      if xi>NBBinsX-1: xi = xi-(NBBinsX-1)           
                      if yi<0: yi = NBBinsY-1+yi
                      if yi>NBBinsX-1: yi = yi-(NBBinsY-1) 
                      # add species if place is empty
                      if model == 'ballisticnormal':
                          if bins[xi,yi] == 0:                         
                              bins[xi,yi] = color
                              hprofiley[xi] = yi
                      if model == 'ballisticcosine':
                          if bins[xi,yi] == 0:                         
                              bins[xi,yi] = color
                              hprofiley[xi] = yi
                      if model == 'surfacediffusion':
                          # find adjacent place within 
                          # surfacediffusion length to
                          # deposit particle
                          refh = hprofiley[xi]
                          xnew = xi
                          for i in range(1,surfacesteps+1,1):
                             xright =  xi+i
                             xleft = xi-i
                             if xright>NBBinsX-1: xright = xright-(NBBinsX-1)
                             if xleft<0: xleft = (NBBinsX-1)-abs(xleft)
                             if hprofiley[xleft]<refh: xnew = xleft            
                             if hprofiley[xright]<refh: xnew = xright
                          hprofiley[xnew] = hprofiley[xnew]+1
                          bins[xnew,int(hprofiley[xnew])] = color
                          
                      break
                   
    print('Step: ',int(k*NBParticlesReplot),' of ', int(NBParticles))        
    # material
    
    if NBParticlesReplot*k == NBParticlesChangeColor*ColorPt:
        if ColorPt/2 == int(ColorPt/2):
          color = 1
        else:
          color = 2  
        ColorPt += 1
    
    x = []
    y = []
    for i in range(NBBinsX):
        for j in range(NBBinsY):
            if bins[i,j] == 1:
               x.append(i+0.5)
               y.append(j+0.5)               
    line1.set_xdata(x)
    line1.set_ydata(y) 
   
    linehprofile.set_xdata(hprofilex)
    linehprofile.set_ydata(hprofiley)
    # calculate roughness
    ht.append(k*NBParticlesReplot)
    hmean.append(np.std(hprofiley))
    linehmean.set_xdata(ht)
    linehmean.set_ydata(hmean)
  
    
    x1 = []
    y1 = []
    for i in range(NBBinsX):
        for j in range(NBBinsY):
            if bins[i,j] == 2:
               x1.append(i+0.5)
               y1.append(j+0.5)               
    line2.set_xdata(x1)
    line2.set_ydata(y1) 
    
    fig.suptitle('Particles: ' + "{:.0f}".format(int(k*NBParticlesReplot)))
    
ColorPt = 1
color = 1    
anim = animation.FuncAnimation(fig,animate,interval=1,frames=int(NBParticles/NBParticlesReplot))
anim.save(animationfile,fps=25,dpi=180)


#x = []
#y = []
#for i in range(NBBinsX):
#    for j in range(NBBinsY):
#            if bins[i,j] == 1:
#               x.append(i+0.5)
#               y.append(j+0.5)  
#                
#ax[0].plot(x,y,color='orange',marker='s',markersize=ms,lw=0)[0]    
