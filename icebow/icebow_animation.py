#
# Icebow
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

NBreflectionsinside = 7
anglemax = 2*np.pi
anglesteps = 1200
  
if NBreflectionsinside == 0:
  animationname = 'icebow_0.mp4' 
  tend = 2.5e-7
  
if NBreflectionsinside == 1:
  animationname = 'icebow_1.mp4' 
  tend = 2.5e-7

if NBreflectionsinside == 2:
  animationname = 'icebow_2.mp4'  
  tend = 2.5e-7

if NBreflectionsinside == 3:
  animationname = 'icebow_3.mp4'  
  tend = 3e-7

if NBreflectionsinside == 7:
  animationname = 'icebow_7.mp4'  
  tend = 6e-7

  
R0 = 5
raysNB = 11

t = 0
dt = 5e-10
c0 = 3e8
steps = int(tend/dt) 


# --------------------------------------------------------------------------
# setup plot
# -------------------------------------------------------------------------
xmax = 15
ymax = 10
fig, ax = plt.subplots(1,1,figsize=(8,5.5))
ax.set(xlim=(-xmax,xmax),ylim=(-ymax,ymax))
ax.axis('off')

ystartmax = R0-0.15
deltay = ystartmax/(raysNB+1)

raysplot = []
for i in range(raysNB):
    if i == 0:
      raysplot.append(ax.plot(0,0,color='r',lw=1,label='n = 1.33',alpha=0.5)[0])
    else:
      raysplot.append(ax.plot(0,0,color='r',lw=1,label ='',alpha=0.5)[0])
for i in range(raysNB):
    if i == 0:
      raysplot.append(ax.plot(0,0,color='g',lw=1,label='n = 1.34',alpha=0.5)[0])
    else:
      raysplot.append(ax.plot(0,0,color='g',lw=1,label='',alpha=0.5)[0])
for i in range(raysNB):
    if i == 0:
       raysplot.append(ax.plot(0,0,color='b',lw=1,label='n = 1.35',alpha=0.5)[0])
    else:
       raysplot.append(ax.plot(0,0,color='b',lw=1,label ='',alpha=0.5)[0])
        
# optical axis
ax.plot([-15,15],[0,0],color='g',label='optical axis')

# raindrop
angle = np.linspace(0,2*np.pi,360)
lensx = R0*np.cos(angle)
lensy = R0*np.sin(angle)


infobox = ''
infobox += 'NB reflections inside: '+ "{:.0f}".format(NBreflectionsinside)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,1,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax.transAxes)   
        
ax.legend(fontsize=8,loc=1)



# --------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------
def animate(k):
    
    global t,raysNB,raycrossing,c0,R0
    
    
    # ---------------------------------------------------------------------------
    # Initialize Rays
    # ---------------------------------------------------------------------------
    raysx = []
    raysy = []
    vx = []
    vy = []
    x = []
    y = []
    xn = []
    yn = []
    rayindex = []
    raycrossing = []

    colorNB = 3
    
    # create hexagon
    icex = []
    icey = []
    R0 = 5
    phase = anglemax/anglesteps*k
    print(k,phase)
    
    px1 = R0*np.cos(phase+np.pi*1/6) 
    py1 = R0*np.sin(phase+np.pi*1/6) 
    px2 = R0*np.cos(phase+np.pi*3/6) 
    py2 = R0*np.sin(phase+np.pi*3/6) 
    px3 = R0*np.cos(phase+np.pi*5/6) 
    py3 = R0*np.sin(phase+np.pi*5/6) 
    px4 = R0*np.cos(phase+np.pi*7/6) 
    py4 = R0*np.sin(phase+np.pi*7/6) 
    px5 = R0*np.cos(phase+np.pi*9/6) 
    py5 = R0*np.sin(phase+np.pi*9/6) 
    px6 = R0*np.cos(phase+np.pi*11/6) 
    py6 = R0*np.sin(phase+np.pi*11/6) 
    px7 = R0*np.cos(phase+np.pi*1/6) 
    py7 = R0*np.sin(phase+np.pi*1/6) 
    
    icex = [px1,px2,px3,px4,px5,px6,px7]
    icey = [py1,py2,py3,py4,py5,py6,py7]
    
    poly = [
       (px1,py1),(px2,py2),(px3,py3),(px4,py4),(px5,py5),(px6,py6),(px7,py7)
       ]
    
    poly_patch = plt.Polygon(poly, color='lightblue', ec='b')
    for p in ax.patches:
        p.set_visible(False)
        p.remove()
    ax.add_patch(poly_patch)    
    
    # red beam
    for i in range(raysNB):

        xstart = -15
        ystart = ystartmax - i * deltay
   
        x.append(xstart)
        y.append(ystart)
   
        xn.append(xstart)
        yn.append(ystart)

        vx.append(c0)
        vy.append(0)
   
        raycrossing.append(0)
   
        raysx.append([xstart])
        raysy.append([ystart])

        rayindex.append(1.33)

    # green beam
    for i in range(raysNB):

        xstart = -15
        ystart = ystartmax - i * deltay
   
        x.append(xstart)
        y.append(ystart)
   
        xn.append(xstart)
        yn.append(ystart)

        vx.append(c0)
        vy.append(0)
   
        raycrossing.append(0)
   
        raysx.append([xstart])
        raysy.append([ystart])

        rayindex.append(1.34)

    # blue beam   
    for i in range(raysNB):

        xstart = -15
        ystart = ystartmax - i * deltay

        x.append(xstart)
        y.append(ystart)
   
        xn.append(xstart)
        yn.append(ystart)

        vx.append(c0)
        vy.append(0)
   
        raycrossing.append(0)
   
        raysx.append([xstart])
        raysy.append([ystart])

        rayindex.append(1.35)
   
    t = 0
    #print('Start simulation')

    
    #########################################
    #  propgate trajectories for all rays  
    #########################################
    for m in range(steps):
     
      t += dt      
      for i in range(colorNB*raysNB):
        xn[i] = x[i] + dt*vx[i]
        yn[i] = y[i] + dt*vy[i]

        #  line ray
        a1 = vy[i]/vx[i]
        x1 = x[i]
        y1 = y[i]
        #print(len(icex))
        # check for entering the raindrop
        for j in range(len(icex)-1):
          # line lens
          a2 = (icey[j+1]-icey[j])/(icex[j+1]-icex[j])
          x2 = icex[j]
          y2 = icey[j]
          x0 = (y2-y1 + a1*x1 - a2*x2)/(a1-a2)
          y0 = a1*(x0-x1)+y1       
          if ( ((x0>=icex[j] and x0<=icex[j+1]) 
              or (x0>=icex[j+1] and x0<=icex[j]))
              and
               ((y0>=icey[j] and y0<=icey[j+1]) 
              or (y0>=icey[j+1] and y0<=icey[j]))
                 and ((x0>=x[i] and x0<=xn[i]) or
                      (x0<=x[i] and x0>=xn[i])
                      )
                 and raycrossing[i] == 0):
                    
             #print('entering',i)
                                       
             c = np.sqrt((x1-x2)**2+(y1-y2)**2)
             a = np.sqrt((x1-x0)**2+(y1-y0)**2)
             b = np.sqrt((x2-x0)**2+(y2-y0)**2)               
             gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))              
           
             raycrossing[i] += 1    
             alpha2 = np.arcsin(1/rayindex[i]*np.sin(np.pi/2-gamma))
                       
             dreh = ((np.pi/2-gamma)-alpha2)
               
             vxn = vx[i]*np.cos(dreh)-vy[i]*np.sin(dreh)
             vyn = vx[i]*np.sin(dreh)+vy[i]*np.cos(dreh)
             vx[i] = vxn/rayindex[i]
             vy[i] = vyn/rayindex[i]
               
             x[i] = x0
             y[i] = y0
             raysx[i].append(x[i])
             raysy[i].append(y[i])  
             xn[i] = x[i] + dt*vx[i]
             yn[i] = y[i] + dt*vy[i]

             break
       
        
        # check for internal reflection inside the raindrop
       # for j in range(len(icex)-1):
          # line lens
       #   a2 = (icey[j+1]-icey[j])/(icex[j+1]-icex[j])
       #   x2 = icex[j]
       #   y2 = icey[j]
       #   x0 = (y2-y1 + a1*x1 - a2*x2)/(a1-a2)
       #   y0 = a1*(x0-x1)+y1
          
          elif ( ((x0>=icex[j] and x0<=icex[j+1]) 
              or (x0>=icex[j+1] and x0<=icex[j]))
             and
               ((y0>=icey[j] and y0<=icey[j+1]) 
              or (y0>=icey[j+1] and y0<=icey[j]))
               and ((x0>=x[i] and x0<=xn[i]) 
                    or (x0<=x[i] and x0>=xn[i]) )
                 and (raycrossing[i]>0 and raycrossing[i] < NBreflectionsinside+1)):        
             
             #print('reflection',i)                                          
                   
             c = np.sqrt((x1-x2)**2+(y1-y2)**2)
             a = np.sqrt((x1-x0)**2+(y1-y0)**2)
             b = np.sqrt((x2-x0)**2+(y2-y0)**2)               
             gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))              
           
             raycrossing[i] += 1    
             dreh = +(2*gamma)
             vxn = vx[i]*np.cos(dreh)-vy[i]*np.sin(dreh)
             vyn = vx[i]*np.sin(dreh)+vy[i]*np.cos(dreh)
             vx[i] = vxn
             vy[i] = vyn
               
             x[i] = x0
             y[i] = y0
             raysx[i].append(x[i])
             raysy[i].append(y[i])  
             xn[i] = x[i] + dt*vx[i]
             yn[i] = y[i] + dt*vy[i]

             break
       
        # check for leaving the raindrop
    #    for j in range(len(icex)-1):
          # line lens
    #      a2 = (icey[j+1]-icey[j])/(icex[j+1]-icex[j])
    #      x2 = icex[j]
    #      y2 = icey[j]
    #      x0 = (y2-y1 + a1*x1 - a2*x2)/(a1-a2)
          elif ( ((x0>=icex[j] and x0<=icex[j+1]) 
              or (x0>=icex[j+1] and x0<=icex[j]))
              and
                 ((y0>=icey[j] and y0<=icey[j+1]) 
              or (y0>=icey[j+1] and y0<=icey[j]))
                 and ((x0>=x[i] and x0<=xn[i])
                      or (x0<=x[i] and x0>=xn[i]))
                 and (raycrossing[i] == NBreflectionsinside+1)):        
   
             #print('leaving',i)  
                                     
             y0 = a1*(x0-x1)+y1
                           
             c = np.sqrt((x1-x2)**2+(y1-y2)**2)
             a = np.sqrt((x1-x0)**2+(y1-y0)**2)
             b = np.sqrt((x2-x0)**2+(y2-y0)**2)               
             gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))              
             raycrossing[i] += 1    
             alpha2 = np.arcsin(rayindex[i]*np.sin(np.pi/2-gamma))
             #print('gamma alpha2', gamma, alpha2)
            
             dreh = -((np.pi/2-gamma)-alpha2)
          
             #print('dreh',dreh)  
          
             vxn = vx[i]*np.cos(dreh)-vy[i]*np.sin(dreh)
             vyn = vx[i]*np.sin(dreh)+vy[i]*np.cos(dreh)
             vx[i] = vxn*rayindex[i]
             vy[i] = vyn*rayindex[i]              
               
             x[i] = x0
             y[i] = y0
             raysx[i].append(x[i])
             raysy[i].append(y[i])  
             xn[i] = x[i] + dt*vx[i]
             yn[i] = y[i] + dt*vy[i]

             break
      
       
       
        x[i] = xn[i] 
        y[i] = yn[i]
        raysx[i].append(x[i])
        raysy[i].append(y[i])         
       
   
       
    for i in range(colorNB*raysNB):    
      raysplot[i].set_xdata(raysx[i])
      raysplot[i].set_ydata(raysy[i]) 
 
    
    
            
    #fig.suptitle('time: '+"{:.2f}".format(t/1e-9)+' ns')
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=anglesteps)
anim.save(animationname,fps=25,dpi=300)