#
# Hertz dipole
# Examples for plasma pyhsics lectures
# contourline algorithm from
# Girwidz, R. V. (2016). 
# Visualizing dipole radiation. 
# European Journal of Physics, 37(6), 065206, pp 1-13.
#
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


model = 'hertzfar'
#model = 'hertznear'

if model == 'hertzfar':
  dt = 1e-10
  tend = 50e-9
  lamb = 3
  animationfile  = 'hertzfar_3dline.mp4'
  xmax = 10
if model == 'hertznear':
  dt = 1e-11
  tend = 7e-9
  lamb = 1
  animationfile  = 'hertznear_3dline.mp4'
  lobesperperiod = 24
  xmax = 1.5
  d0 = 0.01


el = 1.6e-19
eps0 = 8.854e-12
c = 299792458
p0 = 100
f = c/lamb
omega = 2*np.pi*f
period = 1/f


def dipole(t):
    z = p0*np.cos(omega*t) 
    x = 0
    y = 0
    return np.array([x,y,z])

def dipoleprime(t):
    z = -p0*omega*np.sin(omega*t) 
    x = 0
    y = 0
    return np.array([x,y,z])

def dipoleprime2(t):
    z = -p0*omega**2*np.cos(omega*t) 
    x = 0
    y = 0
    return np.array([x,y,z])

def Edipole(r,t):
    rabs = np.linalg.norm(r)
    factor = 1#1/(4*np.pi*eps0)
    tr = t - rabs/c
    term1 = -dipoleprime2(tr)*1/(c**2*rabs)
    term2 = np.dot(dipoleprime2(tr),r)*r/(c**2*rabs**3)
    term3 = 3*np.dot(dipoleprime(tr),r)*r/(c*rabs**4)
    term4 = -dipoleprime(tr)*1/(c*rabs**2)  
    term5 = 3*np.dot(r,dipole(tr))*r/(rabs**5)
    term6 = -dipole(tr)/(rabs**3)
    Efeld = factor*(term1+term2+term3+term4+term5+term6)
    return np.array([Efeld[0],Efeld[1],Efeld[2]])

def Hdipole(r,t):
    xp = r[0]
    yp = r[1]
    r2  = xp*xp + yp*yp
    rabs = np.sqrt( r2 )
    tt  = t - rabs/c
   # p   = dipole(tt)
    pd1 = dipoleprime(tt)[2]
    pd2 = dipoleprime2(tt)[2]
    #print(xp,yp,pd1,pd2)
    sina = xp / rabs
    return pd2 / ( c**2 * rabs ) * sina  +  pd1 / ( c * r2 ) * sina
	 

# --------------------------------------------------------------------------
# setup plot
# -------------------------------------------------------------------------
ymax = xmax
zmax = xmax
fig = plt.figure(figsize=(8,4.5))
ax = plt.axes(projection='3d')
ax.set(xlabel="x (m)",ylabel="y (m)")
ax.set(xlim=(-1,2*xmax),ylim=(-ymax,ymax),zlim=(-zmax,zmax))

x0 = lamb/8
nbefeld = 300
xE = np.linspace(x0,2*xmax,nbefeld)
yE = np.zeros(nbefeld)
zE = np.zeros(nbefeld)

nbquivers = 60
xquiversE = np.linspace(x0,2*xmax,nbquivers)

Eline = ax.plot3D(xE,yE,zE,color='r',lw=0.5,alpha=0.8, label ='E field')[0]

quiverefeld = []
for i in range(nbquivers):
   quiverefeld.append(ax.quiver(xquiversE[i], 0, 0, 0, 0, 1, pivot = 'middle', length = 1, color="r", alpha=0.5))

    
x0 = lamb/16
nbhfeld = 300
xH = np.linspace(x0,2*xmax,nbhfeld)
yH = np.zeros(nbhfeld)
zH = np.zeros(nbhfeld)

nbquivers = 60
xquiversH = np.linspace(x0,2*xmax,nbquivers)

Hline = ax.plot3D(xH,yH,zH,color='g',lw=0.5,alpha=0.8, label ='B field')[0]
ax.legend(fontsize=6,loc=1)

quiverhfeld = []
for i in range(nbquivers):
   quiverhfeld.append(ax.quiver(xquiversH[i], 0, 0, 0, 0, 1, pivot = 'middle', length = 1, color="g", alpha=0.5))


#     quiverleft.append(ax.quiver(0, 0, 0, 1, pivot = 'middle', scale = scalearrows, color="g", alpha=0.5))
          

#infobox = ''
#infobox += 'Wavelength: '+ "{:.1f}".format(lamb)+' (m)\n'
#infobox += 'Frequency: '+ "{:.1e}".format(f)+' (1/s)'
#props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
#ax.text2D(0.02,0.98,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)   

#scalearrows = 20
#lobelw = 0.5

dtopx = [0,0]
dtopy = [0,0]
dtopz = [0,0]
dbottomx = [0,0]
dbottomy = [0,0]
dbottomz = [0,0]
dipoletop = ax.plot3D(dtopx,dtopy,dtopz,color='r',lw=2,alpha=0.8, label ='')[0]           
dipolebottom = ax.plot3D(dbottomx,dbottomy,dbottomz,color='b',lw=2,alpha=0.8, label ='')[0]           
dipoletopmarker = ax.plot3D(0,0,0.5,color='r',markersize=2,lw=0,marker='o',alpha=0.8, label ='')[0]           
dipolebottommarker = ax.plot3D(0,0,-0.5,color='b',markersize=2,lw=0,marker='o',alpha=0.8, label ='')[0]           

# --------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------


nb = 0
t = 0
time = 0
Escale = 10

def animate(k):
    
   global t,time,xE,yE,zE,nbefeld
  
#   x0 = 0.5 
   t += dt 
   time += dt
   # propagate existing lobes
 #  r = np.array([x0,0,0])
 #  Ex, Ey, Ez = Edipole(r,t)/Escale
   
   
   #print(Ex, Ey, Ez)
   
  
   for k in range(nbefeld):
       r = np.array([xE[k],0,0])
       zE[k] = Edipole(r,t)[2]/Escale

   Eline.set_xdata(xE)
   Eline.set_ydata(yE)
   Eline.set_3d_properties(zE)

   for k in range(nbhfeld):
       r = np.array([xH[k],0,0])
       #print(Hdipole(r,t))
       yH[k] = Hdipole(r,t)/Escale
       #print(zH[k])

   Hline.set_xdata(xH)
   Hline.set_ydata(yH)
   Hline.set_3d_properties(zH)
    
   for k in range(nbquivers):    
       quiverefeld[k].remove()
       r = np.array([xquiversE[k],0,0])
       Ez = Edipole(r,t)[2]/Escale
       if Ez>0: 
           direction = 1
           quiverefeld[k] = ax.quiver(xquiversE[k], 0, 0, 0, 0, 1, pivot = 'tail', length = abs(Ez), arrow_length_ratio=1/abs(Ez), color="r", alpha=0.5)
       else:
           direction = -1
           quiverefeld[k] = ax.quiver(xquiversE[k], 0, 0, 0, 0, -1, pivot = 'tail', length = abs(Ez), arrow_length_ratio=1/abs(Ez), color="r", alpha=0.5)

   for k in range(nbquivers):    
       quiverhfeld[k].remove()
       r = np.array([xquiversH[k],0,0])
       Hz = Hdipole(r,t)/Escale
       if Hz>0: 
           direction = 1
           quiverhfeld[k] = ax.quiver(xquiversH[k], 0, 0, 0, 1, 0, pivot = 'tail', length = abs(Hz), arrow_length_ratio=1/abs(Hz), color="g", alpha=0.5)
       else:
           direction = -1
           quiverhfeld[k] = ax.quiver(xquiversH[k], 0, 0, 0, -1, 0, pivot = 'tail', length = abs(Hz), arrow_length_ratio=1/abs(Hz), color="g", alpha=0.5)

   px,py,pz = dipole(t)
   ptop = Escale/p0*pz
   pbottom = -ptop
   dipoletopmarker.set_xdata([0])
   dipoletopmarker.set_ydata([0])
   dipoletopmarker.set_3d_properties([ptop])
   dipolebottommarker.set_xdata([0])
   dipolebottommarker.set_ydata([0])
   dipolebottommarker.set_3d_properties([pbottom])
       
   dtopx = [0,0]
   dtopy = [0,0]
   dtopz = [0,ptop]
   dbottomx = [0,0]
   dbottomy = [0,0]
   dbottomz = [0,pbottom]
   dipoletop.set_xdata(dtopx)
   dipoletop.set_ydata(dtopy)
   dipoletop.set_3d_properties(dtopz)
   dipolebottom.set_xdata(dbottomx)
   dipolebottom.set_ydata(dbottomy)
   dipolebottom.set_3d_properties(dbottomz)

  
       #quiverefeld[k].set_UVC(0,0,0,0,0,1)
 
    
#   px,py,pz = dipole(t)  
#   ptop = 0.5/p0*pz
#   pbottom = -ptop
#   dipoletopmarker.set_xdata([0])
#   dipoletopmarker.set_ydata([ptop])
#   dipolebottommarker.set_xdata([0])
#   dipolebottommarker.set_ydata([pbottom])
       
#   dtopx = [0,0]
#   dtopy = [0,ptop]#
#   dbottomx = [0,0]
#   dbottomy = [0,pbottom]
#   dipoletop.set_xdata(dtopx)
#   dipoletop.set_ydata(dtopy)
#   dipolebottom.set_xdata(dbottomx)
#   dipolebottom.set_ydata(dbottomy)

       
   time = 0
   print("{:.2e}".format(t)+' ns of',' at ',"{:.2e}".format(tend)+' ns')
  
   fig.suptitle('time: '+"{:.1f}".format(t/1e-9)+' ns')

       
  
    
anim = animation.FuncAnimation(fig,animate,interval=1,frames=int(tend/dt))
anim.save(animationfile,fps=25,dpi=300)

