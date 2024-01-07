#
# Rainbow
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

NBreflectionsinside = 7

if NBreflectionsinside == 1:
  animationname = 'rainbow_1.gif' 
  tend = 2.5e-7
if NBreflectionsinside == 2:
  animationname = 'rainbow_2.gif'  
  tend = 2.5e-7
if NBreflectionsinside == 3:
  animationname = 'rainbow_3.gif'  
  tend = 3e-7
if NBreflectionsinside == 7:
  animationname = 'rainbow_7.gif'  
  tend = 6e-7
  
R0 = 5
raysNB = 11

t = 0
dt = 5e-10
c = 3e8
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
      raysplot.append(ax.plot(0,0,color='r',lw=0.6,label ='n$_{water}$ = 1.33',alpha=0.5)[0])
    else:
      raysplot.append(ax.plot(0,0,color='r',lw=0.6,label ='',alpha=0.5)[0])
for i in range(raysNB):
    if i == 0:
      raysplot.append(ax.plot(0,0,color='g',lw=0.6,label ='n$_{water}$ = 1.34',alpha=0.5)[0])
    else:
      raysplot.append(ax.plot(0,0,color='g',lw=0.6,label ='',alpha=0.5)[0])
for i in range(raysNB):
    if i == 0:
       raysplot.append(ax.plot(0,0,color='b',lw=0.6,label ='n$_{water}$ = 1.35',alpha=0.5)[0])
    else:
       raysplot.append(ax.plot(0,0,color='b',lw=0.6,label ='',alpha=0.5)[0])
        
# optical axis
ax.plot([-15,15],[0,0],color='g',label='optical axis')

# raindrop
angle = np.linspace(0,2*np.pi,360)
lensx = R0*np.cos(angle)
lensy = R0*np.sin(angle)
ax.fill_between(lensx,lensy,fc='lightblue')

infobox = ''
infobox += 'NB reflections inside: '+ "{:.0f}".format(NBreflectionsinside)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,1,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax.transAxes)   
        
ax.legend(fontsize=8,loc=1)


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

# red beam
for i in range(raysNB):

   xstart = -15
   ystart = ystartmax - i * deltay
   
   x.append(xstart)
   y.append(ystart)
   
   xn.append(xstart)
   yn.append(ystart)

   vx.append(c)
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

   vx.append(c)
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

   vx.append(c)
   vy.append(0)
   
   raycrossing.append(0)
   
   raysx.append([xstart])
   raysy.append([ystart])

   rayindex.append(1.35)
   
print('Start simulation')

# --------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------
def animate(k):
    
    global t,raysNB,raycrossing
    
    # propgate trajectories for all rays  
    for i in range(colorNB*raysNB):
       xn[i] = x[i] + dt*vx[i]
       yn[i] = y[i] + dt*vy[i]

       # line ray
       a1 = vy[i]/vx[i]
       x1 = x[i]
       y1 = y[i]
      
       # check for entering the raindrop
       if (( xn[i]**2+yn[i]**2<= R0**2) and 
           ( x[i]**2 + y[i]**2>= R0**2) and
           raycrossing[i] == 0):
           
           print('entering',i)
           
           x01 = (a1**2*x1 - a1*y1 - 
                 np.sqrt(R0**2 + a1**2*R0**2 - a1**2*x1**2 + 2*a1*x1*y1 - y1**2))/(1 + a1**2)
           x02 = (a1**2*x1 - a1*y1 + 
                 np.sqrt(R0**2 + a1**2*R0**2 - a1**2*x1**2 + 2*a1*x1*y1 - y1**2))/(1 + a1**2)
           
           if (x01>=x[i] and x01<=xn[i]):
               x0 = x01
           elif x02>=x[i] and x02<=xn[i]: 
               x0 = x02
           #else: break
       
           if ((x01>=x[i] and x01<=xn[i]) or
              (x02>=x[i] and x02<=xn[i])): 
                          
             y0 = a1*(x0-x1)+y1
           
             x2 = x0*1.1
             y2 = y0*1.1
                   
             c = np.sqrt((x1-x2)**2+(y1-y2)**2)
             a = np.sqrt((x1-x0)**2+(y1-y0)**2)
             b = np.sqrt((x2-x0)**2+(y2-y0)**2)               
             gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))              
           
             raycrossing[i] += 1    
             alpha2 = np.arcsin(1/rayindex[i]*np.sin(gamma))
             dreh = -(gamma-alpha2)
               
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
       if (( xn[i]**2 + yn[i]**2>= R0**2) and 
           ( x[i]**2 + y[i]**2<= R0**2) and
           ( raycrossing[i]>0 and raycrossing[i]<NBreflectionsinside+1)):
           print('reflection',i)        
           x01 = (a1**2*x1 - a1*y1 -	 
                 np.sqrt(R0**2 + a1**2*R0**2 - a1**2*x1**2 + 2*a1*x1*y1 - y1**2))/(1 + a1**2)
           x02 = (a1**2*x1 - a1*y1 +	 
                 np.sqrt(R0**2 + a1**2*R0**2 - a1**2*x1**2 + 2*a1*x1*y1 - y1**2))/(1 + a1**2)
           
          
           if (xn[i]<=x01<=x[i] or xn[i]>=x01>=x[i]):
               x0 = x01
           elif (xn[i]<=x02<=x[i] or xn[i]>=x02>=x[i]): 
               x0 = x02
          
           if ((xn[i]<=x01<=x[i] or xn[i]>=x01>=x[i]) or
               (xn[i]<=x02<=x[i] or xn[i]>=x02>=x[i])): 
        
           
        #   if (x01>=x[i] and x01<=xn[i]): 
        #       x0 = x01
        #   elif (x02>=x[i] and x02<=xn[i]):
        #       x0 = x02
          
        #   if ((x01>=x[i] and x01<=xn[i]) or 
        #       (x02>=x[i] and x02<=xn[i])):
                         
             y0 = a1*(x0-x1)+y1
           
             x2 = x0*0.9
             y2 = y0*0.9
                   
             c = np.sqrt((x1-x2)**2+(y1-y2)**2)
             a = np.sqrt((x1-x0)**2+(y1-y0)**2)
             b = np.sqrt((x2-x0)**2+(y2-y0)**2)               
             gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))              
           
             raycrossing[i] += 1    
             dreh = -(np.pi-2*gamma)
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
       if (( xn[i]**2+yn[i]**2>= R0**2) and 
           ( x[i]**2 + y[i]**2<= R0**2) and
           raycrossing[i] == NBreflectionsinside+1):
           print('leaving',i)  
           
           x01 = (a1**2*x1 - a1*y1 - 
                 np.sqrt(R0**2 + a1**2*R0**2 - a1**2*x1**2 + 2*a1*x1*y1 - y1**2))/(1 + a1**2)
           x02 = (a1**2*x1 - a1*y1 + 
                 np.sqrt(R0**2 + a1**2*R0**2 - a1**2*x1**2 + 2*a1*x1*y1 - y1**2))/(1 + a1**2)
           #print('x xn x01 x02',x[i],xn[i],x01,x02)
           if (xn[i]<=x01<=x[i] or xn[i]>=x01>=x[i]):
               x0 = x01
           elif (xn[i]<=x02<=x[i] or xn[i]>=x02>=x[i]): 
               x0 = x02
          
           if ((xn[i]<=x01<=x[i] or xn[i]>=x01>=x[i]) or
               (xn[i]<=x02<=x[i] or xn[i]>=x02>=x[i])): 
            
          #if ((x01>=x[i] and x01<=xn[i]) or
          #     (x02>=x[i] and x02<=xn[i])): 
                          
             y0 = a1*(x0-x1)+y1
           
             x2 = x0*0.9
             y2 = y0*0.9
                   
             c = np.sqrt((x1-x2)**2+(y1-y2)**2)
             a = np.sqrt((x1-x0)**2+(y1-y0)**2)
             b = np.sqrt((x2-x0)**2+(y2-y0)**2)               
             gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))              
             raycrossing[i] += 1    
             alpha2 = np.arcsin(rayindex[i]*np.sin(gamma))
             #print('gamma alpha2', gamma, alpha2)
            
             dreh = (gamma-alpha2)
          
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
       raysplot[i].set_xdata(raysx[i])
       raysplot[i].set_ydata(raysy[i]) 
 
    
    t += dt    
            
    #fig.suptitle('time: '+"{:.2f}".format(t/1e-9)+' ns')
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)