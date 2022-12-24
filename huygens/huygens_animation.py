#
# Reflection thin films
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as colormap

model = "Huygens"

# Initialize Simulation
if model == "Huygens":
    animationname = 'huygens.mp4'
    angledeg = 40 # angle of incidence to normal
    xstart = -14
    xoffset = 1
    raysNB = 10

t = 0
dt = 1e-10
tend = 100e-9
c = 3e8
steps = int(tend/dt) 

# --------------------------------------------------------------------------
# Fresnel equations
# -------------------------------------------------------------------------
def snellius(n1,n2,angle):
    return np.arcsin(n1/n2*np.sin(angle))

def rp(n1,n2,angle1,angle2):
    return (n1*np.cos(angle1)-n2*np.cos(angle2))/(n1*np.cos(angle1)+n2*np.cos(angle2))

def tp(n1,n2,angle1,angle2):
    return (2*n1*np.cos(angle1))/(n2*np.cos(angle1)+n1*np.cos(angle2))

def R(n1,n2,angle1,angle2):
    return rp(n1,n2,angle1,angle2)**2 

def T(n1,n2,angle1,angle2):
    return n2*np.cos(angle2)/(n1*np.cos(angle1))*tp(n1,n2,angle1,angle2)**2 
     

# --------------------------------------------------------------------------
# setup plot
# -------------------------------------------------------------------------
xmax = 15
ymax = 10
fig, ax = plt.subplots(1,1,figsize=(8,4.5))
cmap = colormap.get_cmap('coolwarm')
fig.colorbar(colormap.ScalarMappable(norm=None, cmap=cmap),ax=ax,label='intensity')
raysplot = []
circletop = []
circlebottom = []
#scattermaxima = ax.plot(0,0,color='m',lw=0,markersize=3,label ='maxima',alpha=0.8,marker='o',fillstyle='full')[0]
#scatterminima = ax.plot(0,0,color='g',lw=0,markersize=3,label ='minima',alpha=0.8,marker='o',fillstyle='full',markerfacecolor='white')[0]
#ax.legend(fontsize=6,loc=1)

class layer:
    def __init__(self,y,n,v,angle):
        self.y = y # net charge
        self.n = n # positive charge
        self.v = v # positive charge
        self.angle = angle # negative charge
  
layers = []
NBlayers = 3
for i in range(NBlayers):
  layers.append(layer(y=0,n=1,v=c,angle=0))

layers[0].y = ymax
layers[0].n = 1
layers[0].v = c/layers[0].n
layers[0].angle = angledeg/180*np.pi

layers[1].y = 0
layers[1].n = 2
layers[1].v = c/layers[1].n
layers[1].angle = snellius(layers[0].n,layers[1].n,layers[0].angle)

layers[2].y = -ymax
layers[2].n = 1
layers[2].v = c/layers[2].n
layers[2].angle = snellius(layers[1].n,layers[2].n,layers[1].angle)

for j in range(1,NBlayers-1):
    b3x = [-xmax,xmax]
    b3y = [layers[j].y,layers[j].y]
    b3y2 = [layers[j+1].y,layers[j+1].y]
    if j/2 == int(j/2):
      ax.fill_between(b3x,b3y,b3y2,color='orange',lw=2,alpha=0.8) 
    else:
      ax.fill_between(b3x,b3y,b3y2,color='olive',lw=2,alpha=0.8) 
    
for j in range(1,NBlayers):    
    infobox = ''
    infobox += 'n: '+ "{:.1f}".format(layers[j-1].n)
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
    ax.text(0.02,(layers[j-1].y-(-ymax))/(2*ymax)-0.02,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax.transAxes)   
        
ax.set(xlabel="x (m)",ylabel="y (m)")
ax.set(xlim=(-xmax,xmax),ylim=(-ymax,ymax))



# ---------------------------------------------------------------------------
# Initialize Rays
# ---------------------------------------------------------------------------
x = []
y = []
xn = []
yn = []
phase = []
Int = []
vx = []
vy = []
raysx = []
raysy = []
raysI = []
raysPhase = []

ystart = ymax

for j in range(raysNB):
  x.append(xstart+j*xoffset*np.cos(layers[0].angle))
  y.append(ystart-(raysNB-j)*xoffset*np.sin(layers[0].angle))
  xn.append(xstart+j*xoffset*np.cos(layers[0].angle))
  yn.append(ystart-(raysNB-j)*xoffset*np.sin(layers[0].angle))
  phase.append(0)
  Int.append(1)
  vx.append(layers[0].v*np.sin(layers[0].angle))
  vy.append(-layers[0].v*np.cos(layers[0].angle))
  raysx.append([x[0]])
  raysy.append([y[0]])
  raysI.append(Int[0])
  raysPhase.append(0)
  raycolor = cmap(0.999)
  raysplot.append(ax.plot(0,0,color=raycolor,lw=2,label ='',alpha=0.8)[0])

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
      xc.append(x+r*np.cos(j/90*np.pi))  
      yc.append(y+r*np.sin(j/90*np.pi))
    return xc, yc  

def createcirclebottom(x,y,r):
    xc = []
    yc = []
    for j in range(91):
      xc.append(x+r*np.cos(j/90*np.pi+np.pi))  
      yc.append(y+r*np.sin(j/90*np.pi+np.pi))
    return xc, yc  

# --------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------
def animate(k):
    
    global t,raysNB,timemaxima,maximatrue,circletopNB,circlebottomNB
    
    # propgate trajectories and search for beam split 
    for i in range(raysNB):
       xn[i] = x[i] + dt*vx[i]
       yn[i] = y[i] + dt*vy[i] 
       if x[i]>=-xmax and x[i]<=xmax and y[i]>=-ymax and y[i]<=ymax:
          # search interfaces
          for j in range(NBlayers):
            
              # interface transferred from top to bottom
              if yn[i] < layers[j].y and y[i]>layers[j].y and j==1: 
                # reflect original beam                  
                print('Splitting beam ',i,' at interface: ',j,' at t=',"{:.1f}".format(t/1e-9),' ns of ',"{:.0f}".format(tend/1e-9))
                      
                # create reflected beam
                rayI = raysI[i]*R(layers[j-1].n,layers[j].n,layers[j-1].angle,layers[j].angle) 
                raycolor = cmap(rayI)               
                raysplot.append(ax.plot(0,0,color=raycolor,lw=2,alpha=0.8, label ='')[0])
                raysI.append(rayI)                            
                raysx.append([x[i]])
                raysy.append([y[i]])
                x.append(x[i])
                y.append(y[i])
                xn.append(x[i])
                yn.append(y[i])
                vx.append(layers[j-1].v*np.sin(layers[j-1].angle))
                vy.append(layers[j-1].v*np.cos(layers[j-1].angle))
                raysNB += 1
                
                # create reflected circle
                circletopx.append([x[i]])
                circletopy.append([y[i]])
                circletopr.append(0)
                circletopv.append(layers[j-1].v)
                circletop.append(ax.plot(0,0,color=raycolor,lw=1,alpha=0.8, linestyle = 'dashed', label ='')[0])             
                circletopNB += 1             
                                
                # create transmitted beam                
                rayI = raysI[i]*T(layers[j-1].n,layers[j].n,layers[j-1].angle,layers[j].angle)
                raycolor = cmap(rayI)                               
                raysplot.append(ax.plot(0,0,color=raycolor,lw=2,alpha=0.8, label ='')[0])
                raysI.append(rayI)                            
                raysx.append([xn[i]])
                raysy.append([yn[i]])
                x.append(xn[i])
                y.append(yn[i])
                xn.append(xn[i])
                yn.append(yn[i])
                vx.append(layers[j].v*np.sin(layers[j].angle))
                vy.append(-layers[j].v*np.cos(layers[j].angle))
                raysNB += 1

                # create reflected circle
                circlebottomx.append(xn[i])
                circlebottomy.append(yn[i])
                circlebottomr.append(0)
                circlebottomv.append(layers[j].v)
                circlebottom.append(ax.plot(0,0,color=raycolor,lw=1,alpha=0.8, linestyle='dashed',label ='')[0])             
                circlebottomNB += 1             
                
                # stop propagation original beam
                vy[i] = 0 
                vx[i] = 0           
                     
    # update rays     
    for i in range(raysNB):   
            x[i] = xn[i] 
            y[i] = yn[i]
            raysx[i].append(x[i])
            raysy[i].append(y[i])
            raysplot[i].set_xdata(raysx[i])
            raysplot[i].set_ydata(raysy[i]) 
            
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
    fig.suptitle('time: '+"{:.2f}".format(t/1e-9)+' ns')
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)