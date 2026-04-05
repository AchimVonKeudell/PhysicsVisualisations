#
# Mirror Optics
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

animationname = 'mirror.mp4'  
lensx = []
lensy = []
lensx2 = []
lensy2 = []

lensz = [5]
lensR = [3]
nlens = [3]
lensangle = [np.pi/3]
nlensdirection = [1]
nbsphere = 240
      
t = 0
timemaxima = 0 
dt = 1e-10
tend = 1.7e-7
c = 3e8
steps = int(tend/dt) 

  
# --------------------------------------------------------------------------
# setup plot
# -------------------------------------------------------------------------
xmax = 15
ymax = 10
fig, ax = plt.subplots(1,1,figsize=(8,4.5))

raysNB = 17
ystartmax = 7
deltay = 2*ystartmax/(raysNB-1)

raysplot = []
for i in range(raysNB):
  raysplot.append(ax.plot(0,0,color='b',lw=0.7,label ='',alpha=0.8)[0])

raysplot2 = []
for i in range(raysNB):
  raysplot2.append(ax.plot(0,0,color='r',lw=0.7,label ='',linestyle='dashed',alpha=0.8)[0])

# optical axis
ax.plot([-15,15],[0,0],color='g',label='optical axis')

# mirror 
x = np.linspace(-6,5,nbsphere)
lensspherex = x
lensspherey = 0+lensR[0]*np.sqrt(abs(x-5))
ax.fill_between(lensspherex,lensspherey,ymax,fc='lightblue',ec='b')
x = np.linspace(5,-6,nbsphere)
lensspherex2 = x
lensspherey2 = 0-lensR[0]*np.sqrt(abs(x-5))
ax.fill_between(lensspherex2,lensspherey2,-ymax,fc='lightblue',ec='b')
lensx.append([*lensspherex,*lensspherex2])
lensy.append([*lensspherey,*lensspherey2])

ax.fill_between([5,xmax],[ymax,ymax],color='lightblue')    
ax.fill_between([5,xmax],[-ymax,-ymax],color='lightblue')    

# circle 
angle = np.linspace(-np.pi/2,np.pi/2,nbsphere)
lensspherex = -5+10*np.cos(angle)
lensspherey = 10*np.sin(angle)
ax.plot(lensspherex,lensspherey,lw=1,color='r',linestyle='dashed')
lensx.append(lensspherex)
lensy.append(lensspherey)


#infobox = ''
#infobox += 'n$_{medium}$: '+ "{:.1f}".format(nlens[0])
#props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
#ax.text(0.02,1,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax.transAxes)   
        
#ax.set(xlabel="x (m)",ylabel="y (m)")
ax.set(xlim=(-xmax,xmax),ylim=(-ymax,ymax))
ax.axis('off')
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

raysx2 = []
raysy2 = []
vx2 = []
vy2 = []
xp2 = []
yp2 = []
xn2 = []
yn2 = []

for i in range(raysNB):

   xstart = -15
   ystart = ystartmax - i * deltay
   
   if ystart == 0: ystart = ystartmax

   x.append(xstart)
   y.append(ystart)
   
   xn.append(xstart)
   yn.append(ystart)

   vx.append(c)
   vy.append(0)
   
   raysx.append([x[i]])
   raysy.append([y[i]])
   
for i in range(raysNB):

   xstart = -15
   ystart = ystartmax - i * deltay
   
   if ystart == 0: ystart = ystartmax

   xp2.append(xstart)
   yp2.append(ystart)
   
   xn2.append(xstart)
   yn2.append(ystart)

   vx2.append(c)
   vy2.append(0)
   
   raysx2.append([xstart])
   raysy2.append([ystart])   


print('Start simulation')
maximatrue = True

# --------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------
def animate(k):
    
    global t,raysNB,timemaxima,maximatrue
    
    # propgate trajectories for rays mirror
    for i in range(raysNB):
       xn[i] = x[i] + dt*vx[i]
       yn[i] = y[i] + dt*vy[i]

       # line ray
       a1 = vy[i]/vx[i]
       x1 = x[i]
       y1 = y[i]
       
       # check for crossing of mirror
       for j in range(len(lensx[0])-1):
           # line lens
           a2 = (lensy[0][j+1]-lensy[0][j])/(lensx[0][j+1]-lensx[0][j])
           x2 = lensx[0][j]
           y2 = lensy[0][j]
           x0 = (y2-y1 + a1*x1 - a2*x2)/(a1-a2)
           if ( ((lensx[0][j]<=x0<=lensx[0][j+1]) 
              or (lensx[0][j+1]<=x0<=lensx[0][j]))
                 and (x[i]<=x0<=xn[i] or x[i]>=x0>=xn[i])):
               
               y0 = a1*(x0-x1)+y1
               c = np.sqrt((x1-x2)**2+(y1-y2)**2)
               a = np.sqrt((x1-x0)**2+(y1-y0)**2)
               b = np.sqrt((x2-x0)**2+(y2-y0)**2)
               #print(a,b,c)
               gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))
               dreh = -((np.pi-2*(np.pi/2-gamma)))
               
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

               #print('crossing at ', x0,x[i],xn[i])
               #print('v ', vx[i],vy[i])
               break 

       x[i] = xn[i] 
       y[i] = yn[i]
       raysx[i].append(x[i])
       raysy[i].append(y[i])         
       raysplot[i].set_xdata(raysx[i])
       raysplot[i].set_ydata(raysy[i]) 
       
    # propgate trajectories for rays sphere
    for i in range(raysNB):
       xn2[i] = xp2[i] + dt*vx2[i]
       yn2[i] = yp2[i] + dt*vy2[i]

       # line ray
       a1 = vy2[i]/vx2[i]
       x1 = xp2[i]
       y1 = yp2[i]
       
       # check for crossing of mirror
       for j in range(len(lensx[1])-1):
           # line lens
           a2 = (lensy[1][j+1]-lensy[1][j])/(lensx[1][j+1]-lensx[1][j])
           x2 = lensx[1][j]
           y2 = lensy[1][j]
           x0 = (y2-y1 + a1*x1 - a2*x2)/(a1-a2)
           if ( ((lensx[1][j]<=x0<=lensx[1][j+1]) 
              or (lensx[1][j+1]<=x0<=lensx[1][j]))
                 and (xp2[i]<=x0<=xn2[i] or xp2[i]>=x0>=xn2[i])):
               
               y0 = a1*(x0-x1)+y1
               c = np.sqrt((x1-x2)**2+(y1-y2)**2)
               a = np.sqrt((x1-x0)**2+(y1-y0)**2)
               b = np.sqrt((x2-x0)**2+(y2-y0)**2)
               #print(a,b,c)
               gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))
               dreh = ((np.pi-2*(np.pi/2-gamma)))
               
               vxn = vx2[i]*np.cos(dreh)-vy2[i]*np.sin(dreh)
               vyn = vx2[i]*np.sin(dreh)+vy2[i]*np.cos(dreh)
               vx2[i] = vxn
               vy2[i] = vyn
               
               xp2[i] = x0
               yp2[i] = y0
               raysx2[i].append(xp2[i])
               raysy2[i].append(yp2[i])  
               xn2[i] = xp2[i] + dt*vx2[i]
               yn2[i] = yp2[i] + dt*vy2[i]

               #print('crossing at ', x0,x[i],xn[i])
               #print('v ', vx[i],vy[i])
               break 

       xp2[i] = xn2[i] 
       yp2[i] = yn2[i]
       raysx2[i].append(xp2[i])
       raysy2[i].append(yp2[i])         
       raysplot2[i].set_xdata(raysx2[i])
       raysplot2[i].set_ydata(raysy2[i]) 

    
    t += dt    
    
            
    #fig.suptitle('time: '+"{:.2f}".format(t/1e-9)+' ns')
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)