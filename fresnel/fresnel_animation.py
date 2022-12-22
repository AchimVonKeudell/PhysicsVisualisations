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

model = "Airy"
#model = "Langevin"

# Initialize Simulation
if model == "Airy":
    animationname = 'fresnel_airy.mp4'

t = 0
timemaxima = 0 
dt = 1e-10
dtplot = 2e-10
dtmaxima = 3e-9
tend = 800e-9
c = 3e8
steps = int(tend/dtplot) 

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
xmax = 18
ymax = 10
fig, ax = plt.subplots(1,1,figsize=(8,4.5))
cmap = colormap.get_cmap('coolwarm')
fig.colorbar(colormap.ScalarMappable(norm=None, cmap=cmap),ax=ax,label='intensity')
raysplot = []
raycolor = cmap(0.999)
raysplot.append(ax.plot(0,0,color=raycolor,lw=2,label ='',alpha=0.8)[0])
#scattermaxima = ax.plot(0,0,color='m',lw=0,label ='',alpha=0.8,marker='o')[0]

class layer:
    def __init__(self,y,n,v,angle):
        self.y = y # net charge
        self.n = n # positive charge
        self.v = v # positive charge
        self.angle = angle # negative charge
  
layers = []
NBlayers = 4
for i in range(NBlayers):
  layers.append(layer(y=0,n=1,v=c,angle=0))

layers[0].y = ymax
layers[0].n = 1
layers[0].v = c/layers[0].n
layers[0].angle = 20/180*np.pi

layers[1].y = 2
layers[1].n = 2
layers[1].c = c/layers[1].n
layers[1].angle = snellius(layers[0].n,layers[1].n,layers[0].angle)

layers[2].y = -2
layers[2].n = 3
layers[2].v = c/layers[2].n
layers[2].angle = snellius(layers[1].n,layers[2].n,layers[1].angle)

layers[3].y = -ymax
layers[3].n = 1
layers[3].v = c/layers[3].n
layers[3].angle = snellius(layers[2].n,layers[3].n,layers[2].angle)

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
raysx = []
raysy = []
raysI = []
raysNB = 1

x = [-12]
y = [ymax]
xn = [-12]
yn = [ymax]
Int = [1]
vx = [layers[0].v*np.sin(layers[0].angle)]
vy = [-layers[0].v*np.cos(layers[0].angle)]
raysx.append([x[0]])
raysy.append([y[0]])
raysI.append(Int[0])

#maximax = [-12]
#maximay = [ymax]
#maximax.append([-12])
#maximay.append([ymax])


print('Start simulation')

# --------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------
def animate(k):
    
    global t,raysNB,timemaxima
    
    # propgate trajectories for all rays within dtreplot
    
    for i in range(raysNB):
       time = 0       
       beamsplitting = False
       # propgate each rays for dtreplot          
       while x[i]>=-xmax and x[i]<=xmax and y[i]>=-ymax and y[i]<=ymax and time<dtplot:
          time += dt      
          # propagate rays        
          xn[i] = x[i] + dt*vx[i]
          yn[i] = y[i] + dt*vy[i]
          # search interfaces
          for j in range(NBlayers):
            
              # interface transferred from top to bottom
              if yn[i] < layers[j].y and y[i]>layers[j].y: 
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
                
                # stop propagation original beam
                vy[i] = 0 
                vx[i] = 0           
           
                beamsplitting = True         
         
              # interface transferred from bottom to top
              if yn[i] > layers[j].y and y[i]<layers[j].y and layers[j].y != ymax: 
                # reflect original beam                  
                print('Splitting beam ',i,' at interface: ',j,' at t=',"{:.2f}".format(t/1e-9),' ns of ',"{:.0f}".format(tend/1e-9))
                
                rayI = raysI[i]*R(layers[j].n,layers[j-1].n,layers[j].angle,layers[j-1].angle)              
                raycolor = cmap(rayI)                               
                raysplot.append(ax.plot(0,0,color=raycolor,lw=2,alpha=0.8, label ='')[0])
                raysI.append(rayI)                            
                raysx.append([x[i]])
                raysy.append([y[i]])
                x.append(x[i])
                y.append(y[i])
                xn.append(x[i])
                yn.append(y[i])
                vx.append(layers[j].v*np.sin(layers[j].angle))
                vy.append(-layers[j].v*np.cos(layers[j].angle))
                raysNB += 1
                                
                # create transmitted beam                                
                rayI = raysI[i]*T(layers[j].n,layers[j-1].n,layers[j].angle,layers[j-1].angle)              
                raycolor = cmap(rayI)                               
                
                raysplot.append(ax.plot(0,0,color=raycolor,lw=2,alpha=0.8, label ='')[0])
                raysI.append(rayI)                            
                raysx.append([xn[i]])
                raysy.append([yn[i]])
                x.append(xn[i])
                y.append(yn[i])
                xn.append(xn[i])
                yn.append(yn[i])
                vx.append(layers[j-1].v*np.sin(layers[j-1].angle))
                vy.append(layers[j-1].v*np.cos(layers[j-1].angle))
                raysNB += 1
                
                # stop propagation original beam
                vy[i] = 0 
                vx[i] = 0           
           
                beamsplitting = True
            
         
          if beamsplitting == False:   
            x[i] = xn[i] 
            y[i] = yn[i]
            raysx[i].append(x[i])
            raysy[i].append(y[i])
          
            raysplot[i].set_xdata(raysx[i])
            raysplot[i].set_ydata(raysy[i]) 
          else:
            break 
    # update time after all rays    
    t += time
    timemaxima += time
          
    # create maxima
    #if timemaxima > dtmaxima:
    #   for i in range(raysNB):
    #       maximax.append(x[i])
    #       maximay.append(y[i])
    #   scattermaxima.set_xdata(maximax)    
    #   scattermaxima.set_ydata(maximay)    
    #   timemaxima = 0
    fig.suptitle('time: '+"{:.2f}".format(t/1e-9)+' ns')
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=250,dpi=300)