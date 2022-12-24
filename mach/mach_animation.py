#
# Mach cones
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as colormap

model = "Mach"

# Initialize Simulation
if model == "Mach":
    animationname = 'mach.mp4'
    m1 = 0.6 # angle of incidence to normal
    m2 = 1.5
    
t = 0
dt = 1e-4
tend = 100e-3
timemaxima = 3e-3
c = 3e2
vtop = m1*c
vbottom = m2*c
steps = int(tend/dt) 

# --------------------------------------------------------------------------
# setup plot
# -------------------------------------------------------------------------
xmax = 18
ymax = 10
fig, ax = plt.subplots(1,1,figsize=(8,4.5))

circletop = []
circlebottom = []
top = ax.plot(0,0,color='g',lw=0,markersize=5,label ='maxima',alpha=0.8,marker='o',fillstyle='full')[0]
bottom = ax.plot(0,0,color='orange',lw=0,markersize=5,label ='minima',alpha=0.8,marker='o',fillstyle='full')[0]
#ax.legend(fontsize=6,loc=1)

    
infobox = ''
infobox += 'Mach number: '+ "{:.1f}".format(m1)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.95,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax.transAxes)   

infobox = ''
infobox += 'Mach number: '+ "{:.1f}".format(m2)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.5,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax.transAxes)   

infobox = ''
infobox += 'sound velocity: '+ "{:.0f}".format(c)+' m/s'
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.98,0.95,infobox, fontsize=8,bbox=props,verticalalignment='top',horizontalalignment='right',transform=ax.transAxes)   

        
ax.set(xlabel="x (m)",ylabel="y (m)")
ax.set(xlim=(-xmax,xmax),ylim=(-ymax,ymax))



# ---------------------------------------------------------------------------
# Initialize Rays
# ---------------------------------------------------------------------------

xtop = -xmax
ytop = 0.5*ymax

xbottom = -xmax
ybottom = -0.5*ymax

circletopx = []
circletopy = []
circletopr = []
circletopv = []
circletopNB = 0

circlebottomx = []
circlebottomy = []
circlebottomr = []
circlebottomv = []
circlebottomNB = 0

print('Start simulation')
maximatrue = True

def createcircletop(x,y,r):
    xc = []
    yc = []
    for j in range(91):
      xc.append(x+r*np.cos(j/90*2*np.pi))  
      yc.append(y+r*np.sin(j/90*2*np.pi))
    return xc, yc  

def createcirclebottom(x,y,r):
    xc = []
    yc = []
    for j in range(91):
      xc.append(x+r*np.cos(j/90*2*np.pi))  
      yc.append(y+r*np.sin(j/90*2*np.pi))
    return xc, yc  

# --------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------
time = 0
t = 0
def animate(k):
    
    global t,time,timemaxima,circletopNB,circlebottomNB
    global vtop,vbottom,xbottom,xtop,ytop,ybottom
    
    # propgate trajectories 
    t += dt
    time += dt
    
    xtop += dt*vtop
    xbottom += dt*vbottom
                     
    # update rays     
    top.set_xdata(xtop)
    top.set_ydata(ytop) 
    
    bottom.set_xdata(xbottom)
    bottom.set_ydata(ybottom) 
    
    if time>timemaxima:
      circletopx.append(xtop)
      circletopy.append(ytop)
      circletopr.append(0)
      circletopv.append(c)
      circletop.append(ax.plot(0,0,color='g',lw=1,alpha=0.8, linestyle = 'dashed', label ='')[0])             
      circletopNB += 1    
    
      circlebottomx.append(xbottom)
      circlebottomy.append(ybottom)
      circlebottomr.append(0)
      circlebottomv.append(c)
      circlebottom.append(ax.plot(0,0,color='orange',lw=1,alpha=0.8, linestyle = 'dashed', label ='')[0])             
      circlebottomNB += 1    
    
      time = 0
            
    # update circles        
    for k in range(circlebottomNB):
              circlebottomr[k] = circlebottomr[k] + dt*circlebottomv[k]
              xc,yc = createcirclebottom(circlebottomx[k],circlebottomy[k],circlebottomr[k])
              circlebottom[k].set_xdata(xc)    
              circlebottom[k].set_ydata(yc)
    for k in range(circletopNB):
              circletopr[k] = circletopr[k] + dt*circletopv[k]
              xc,yc = createcircletop(circletopx[k],circletopy[k],circletopr[k])
              circletop[k].set_xdata(xc)    
              circletop[k].set_ydata(yc)   
     
    t += dt                    
    fig.suptitle('time: '+"{:.2f}".format(t/1e-3)+' ms')
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)