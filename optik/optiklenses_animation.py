#
# Geomertic Optics
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

model = "HalfSphere"
model = "LensFocusing"
model = "LensDeFocusing"
model = "Kepler"
model = "Galilei"

if model == "HalfSphere":
    animationname = 'HalfSphere.mp4'  
    lensx = []
    lensy = []
    lensNB = 1
    lensz = [-6]
    lensR = [13]
    nlens = [3]
    angleincidence = 0
    lensangle = [np.pi/3]
    nlensdirection = ['r']
    nbsphere = 240
    tend = 2.5e-7
    dt = 1e-10
    treverse = tend - dt
    treverseplot = tend - dt
    directionreverse = -1
    reversebeam = False
 
if model == "LensFocusing":
    animationname = 'lens_focusing.mp4'
    lensx = []
    lensy = []
    lensNB = 2
    lensz = [-6,0]
    lensR = [12,-12]
    nlens = [1.5,1.5]
    angleincidence = 0
    lensangle = [np.pi/4*0.92,np.pi/4*0.92]
    nlensdirection = ['r','l']
    nbsphere = 190
    tend = 1.5e-7
    dt = 1e-10
    treverse = 1.3e-7
    treverseplot = 1.4e-7
    directionreverse = -1
    reversebeam = True

if model == "LensDeFocusing":
    animationname = 'lens_defocusing.mp4'
    lensx = []
    lensy = []
    lensNB = 2
    lensz = [-3,0]
    lensR = [-6,6]
    nlens = [1.5,1.5]
    angleincidence = 0
    lensangle = [np.pi/4*0.92,np.pi/4*0.92]
    nlensdirection = ['r','l']
    nbsphere = 190
    tend = 1.5e-7
    dt = 1e-10
    treverse = 1.3e-7
    treverseplot = 1.4e-7
    directionreverse = -1
    reversebeam = True

if model == "Kepler":
    animationname = 'kepler_telescope.mp4'
    lensx = []
    lensy = []
    lensNB = 4
    position2 = 6.5
    lensz = [-13,-7,position2,position2+3.2]
    lensR = [12,-12,4,-4]
    nlens = [1.5,1.5,1.5,1.5]
    angleincidence = 2*np.pi/180
    lensangle = [np.pi/4*0.92,np.pi/4*0.92,np.pi/3*0.91,np.pi/3*0.91]
    nlensdirection = ['r','l','r','l']
    nbsphere = 360
    tend = 1.3e-7
    dt = 2e-10
    treverse = 1.1e-7
    treverseplot = 1.1e-7
    directionreverse = -1
    reversebeam = True

if model == "Galilei":
    animationname = 'galilei_telescope.mp4'
    lensx = []
    lensy = []
    lensNB = 4
    position2 = -3
    lensz = [-13,-7,position2,position2+1]
    lensR = [12,-12,-6,6]
    nlens = [1.5,1.5,1.5,1.5]
    angleincidence = 6*np.pi/180
    lensangle = [np.pi/4*0.92,np.pi/4*0.92,np.pi/3*0.92,np.pi/3*0.92]
    nlensdirection = ['r','l','r','l']
    nbsphere = 360
    tend = 1.3e-7
    dt = 2e-10
    treverse = 0.4e-7
    treverseplot = 1.1e-7
    directionreverse = 1
    reversebeam = True
     

t = 0
timemaxima = 0 


c = 3e8
steps = int(tend/dt) 

  

# --------------------------------------------------------------------------
# setup plot
# -------------------------------------------------------------------------
xmax = 15
ymax = 10
fig, ax = plt.subplots(1,1,figsize=(8,4.5))

raysNB = 9
ystartmax = 1.6
deltay = 2*ystartmax/(raysNB-1)

raysplot = []
for i in range(raysNB):
  raysplot.append(ax.plot(0,0,color='b',lw=1,label ='',alpha=0.8)[0])

raysplot2 = []
for i in range(raysNB):
  raysplot2.append(ax.plot(0,0,color='r',lw=1,label ='',linestyle='dashed',alpha=0.8)[0])

# optical axis
ax.plot([-15,15],[0,0],color='g',label='optical axis')

# lens sphere
for i in range(lensNB):
  if lensR[i]>0: # curvature to the right
    angle = np.linspace(-lensangle[i],0,nbsphere)
    lensspherex = (lensz[i]+lensR[i])-lensR[i]*np.cos(angle)
    lensspherey = 0+lensR[i]*np.sin(angle)
    angle = np.linspace(0,lensangle[i],nbsphere)
    lensspherex2 = (lensz[i]+lensR[i])-lensR[i]*np.cos(angle)
    lensspherey2 = 0+lensR[i]*np.sin(angle)
  else: # curvature to the left
    angle = np.linspace(-lensangle[i],0,nbsphere)
    lensspherex = (lensz[i]-abs(lensR[i]))+abs(lensR[i])*np.cos(angle)
    lensspherey = 0+abs(lensR[i])*np.sin(angle)     
    angle = np.linspace(0,lensangle[i],nbsphere)
    lensspherex2 = (lensz[i]-abs(lensR[i]))+abs(lensR[i])*np.cos(angle)
    lensspherey2 = 0+abs(lensR[i])*np.sin(angle)     
  lensx.append([*lensspherex,*lensspherex2])
  lensy.append([*lensspherey,*lensspherey2])
  
  if nlensdirection[i] == 'r' and lensR[i]>0: # left half of a concave lens
    ax.fill_between(lensspherex,lensspherey,fc='lightblue')
    ax.fill_between(lensspherex2,lensspherey2,fc='lightblue')
    ax.text(lensz[i],-ymax+2,'R='+"{:.0f}".format(lensR[i]), fontsize=8,horizontalalignment='center')   
  if nlensdirection[i] == 'l' and lensR[i]<0: # right half of a concave lens
    ax.fill_between(lensspherex,lensspherey,fc='lightblue')
    ax.fill_between(lensspherex2,lensspherey2,fc='lightblue')
    ax.text(lensz[i],-ymax+2,'R='+"{:.0f}".format(lensR[i]), fontsize=8,horizontalalignment='center')   
  if nlensdirection[i] == 'l' and lensR[i]>0: # left half of a concave lens
    ax.fill_between(lensspherex,lensspherey,-abs(lensR[i]),fc='lightblue')
    ax.fill_between(lensspherex2,lensspherey2,abs(lensR[i]),fc='lightblue')
    ax.text(lensz[i]-abs(0.5*lensR[i]),-ymax+2,'R='+"{:.0f}".format(lensR[i]), fontsize=8,horizontalalignment='center')   
  if nlensdirection[i] == 'r' and lensR[i]<0: # right half of a concave lens
    ax.fill_between(lensspherex,lensspherey,-abs(lensR[i]),fc='lightblue')
    ax.fill_between(lensspherex2,lensspherey2,abs(lensR[i]),fc='lightblue')
    ax.text(lensz[i]+abs(0.5*lensR[i]),-ymax+2,'R='+"{:.0f}".format(lensR[i]), fontsize=8,horizontalalignment='center')   


if model=='HalfSphere':
  ax.fill_between([lensz[0]+0.5*lensR[0],xmax],[ymax,ymax],color='lightblue')    
  ax.fill_between([lensz[0]+0.5*lensR[0],xmax],[-ymax,-ymax],color='lightblue')    

if model=='Galilei':
  ax.fill_between([lensz[2],lensz[3]],[abs(lensR[3]),abs(lensR[3])],color='lightblue')    
  ax.fill_between([lensz[2],lensz[3]],[-abs(lensR[3]),-abs(lensR[3])],color='lightblue')    

if model=='LensDeFocusing':
  ax.fill_between([lensz[0],lensz[1]],[abs(lensR[0]),abs(lensR[0])],color='lightblue')    
  ax.fill_between([lensz[0],lensz[1]],[-abs(lensR[0]),-abs(lensR[0])],color='lightblue')    

infobox = ''
infobox += 'n$_{medium}$: '+ "{:.1f}".format(nlens[0])
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,1,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax.transAxes)   
        
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
vxreverse = []
vyreverse = []
xreverse = []
yreverse = []
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

   xstartn = -15
   ystartn = ystartmax - i * deltay
   
   if (model != 'Kepler' and model != 'Galilei') and ystartn == 0: ystartn = ystartmax
   
   xstart = xstartn*np.cos(angleincidence)-ystartn*np.sin(angleincidence)
   ystart = xstartn*np.sin(angleincidence)+ystartn*np.cos(angleincidence)

   x.append(xstart)
   y.append(ystart)
   xreverse.append(xstart)
   yreverse.append(ystart)
   
   xn.append(xstart)
   yn.append(ystart)

   vx.append(c*np.cos(angleincidence))
   vy.append(c*np.sin(angleincidence))
   vxreverse.append(c*np.cos(angleincidence))
   vyreverse.append(c*np.sin(angleincidence))
   
   raysx.append([x[i]])
   raysy.append([y[i]])


print('Start simulation')
maximatrue = True

# --------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------
def animate(k):
    
    global t,raysNB,timemaxima,maximatrue
    
    # propgate trajectories for all rays  
    for i in range(raysNB):
       xn[i] = x[i] + dt*vx[i]
       yn[i] = y[i] + dt*vy[i]

       # line ray
       a1 = vy[i]/vx[i]
       x1 = x[i]
       y1 = y[i]
       
       if t>treverse and t<=treverse+dt and reversebeam:
         xreverse[i] = x[i]
         yreverse[i] = y[i]
         vxreverse[i] = vx[i]
         vyreverse[i] = vy[i]      
           
       
       # check for crossing of lens
       for k in range(lensNB):
         for j in range(len(lensx[k])-1):
           # line lens
           a2 = (lensy[k][j+1]-lensy[k][j])/(lensx[k][j+1]-lensx[k][j])
           x2 = lensx[k][j]
           y2 = lensy[k][j]
           x0 = (y2-y1 + a1*x1 - a2*x2)/(a1-a2)
           if ( ((x0>=lensx[k][j] and x0<=lensx[k][j+1]) 
              or (x0>=lensx[k][j+1] and x0<=lensx[k][j]))
                 and (x0>=x[i] and x0<=xn[i])):
               
               y0 = a1*(x0-x1)+y1
               c = np.sqrt((x1-x2)**2+(y1-y2)**2)
               a = np.sqrt((x1-x0)**2+(y1-y0)**2)
               b = np.sqrt((x2-x0)**2+(y2-y0)**2)
               #print(a,b,c)
               gamma = np.arccos((a**2+b**2-c**2)/(2*a*b))
               #print('einfall',(np.pi/2-gamma)*180/np.pi)
               if nlensdirection[k]=='r':
                 alpha2 = np.arcsin(1/nlens[k]*np.sin(np.pi/2-gamma))
               else:
                 alpha2 = np.arcsin(nlens[k]*np.sin(np.pi/2-gamma))
               #print('ausfall',alpha2*180/np.pi)
               dreh = -((np.pi/2-gamma)-alpha2)
               #print('dreh',i,j,dreh*180/np.pi)
               
               vxn = vx[i]*np.cos(dreh)-vy[i]*np.sin(dreh)
               vyn = vx[i]*np.sin(dreh)+vy[i]*np.cos(dreh)
                      
               if nlensdirection[k]=='r':             
                 vx[i] = vxn/nlens[k]
                 vy[i] = vyn/nlens[k]
               else:
                 vx[i] = vxn*nlens[k]
                 vy[i] = vyn*nlens[k]                   
               
               x[i] = x0
               y[i] = y0
               raysx[i].append(x[i])
               raysy[i].append(y[i])  
               xn[i] = x[i] + dt*vx[i]
               yn[i] = y[i] + dt*vy[i]

               #print('crossing at ', x0,x[i],xn[i])
               #print('v ', vx[i],vy[i])
               break

       #if x[i]>14: print(np.arctan(vy[i]/vx[i])*180/np.pi)
       x[i] = xn[i] 
       y[i] = yn[i]
       raysx[i].append(x[i])
       raysy[i].append(y[i])         
       raysplot[i].set_xdata(raysx[i])
       raysplot[i].set_ydata(raysy[i]) 
 
    if t>treverseplot and reversebeam:
     for i in range(raysNB):    
       raysplot2[i].set_xdata([xreverse[i],xreverse[i]+directionreverse*vxreverse[i]*tend])
       raysplot2[i].set_ydata([yreverse[i],yreverse[i]+directionreverse*vyreverse[i]*tend]) 
 
    t += dt    
    timemaxima += dt
            
    #fig.suptitle('time: '+"{:.2f}".format(t/1e-9)+' ns')
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)