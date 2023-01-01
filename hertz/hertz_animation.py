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
model = 'hertznear'

if model == 'hertzfar':
  dt = 1e-11
  tend = 20e-9
  lamb = 1
  animationfile  = 'hertzfar.mp4'
  lobesperperiod = 16
  xmax = 4
  d0 = 0.01
if model == 'hertznear':
  dt = 1e-11
  tend = 10e-9
  lamb = 1
  animationfile  = 'hertznear.mp4'
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
dtnewlobe = period/lobesperperiod
dt = period/128

# plotting parameter
mindistance = 0.04 # do not show lobes close to origin
maxstep = 0.2
nbpoynting = 200 # grid field
poyntingscale = 1e5 # scale field to keep values [-1,1]

# ------------------------------------------------------------------------
# electrical dipole
# ------------------------------------------------------------------------
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
    rabs = np.linalg.norm(r)
    factor = 1#1/(4*np.pi*eps0)
    tr = t - rabs/c
    term1 = np.cross(dipoleprime2(tr),r)/(c*rabs**2)
    term2 = np.cross(dipoleprime(tr),r)/(rabs**3)
    Hfeld = factor*(term1+term2)
    return np.array([Hfeld[0],Hfeld[1],Hfeld[2]])

# --------------------------------------------------------------------------
# Decide whether field line is dranw 
# in clockwise or anticlockwise orienttation
# --------------------------------------------------------------------------
def orientationfieldline(t):
    nperiods = np.trunc(t/(period))
    tt = t - nperiods*period
    # field line clockwise
    if (tt > 0) and (tt < period/4 ):
        return 0
    if (tt > period/4 ) and (tt < period/2 ):
        return 1
    # field line counter clockwise
    if  (tt > period/2) and (tt < 3*period/4 ):
        return 0
    if  (tt > 3*period/4) and (tt < period):
        return -1


# --------------------------------------------------------------------------
# Feldlinien
# --------------------------------------------------------------------------
def Efieldline2(rstart,t):

    def rnext(r,t,rechdir):
      E1 = Edipole(r,t)
      Eabs = np.linalg.norm(E1)
      d = d0/np.sqrt(E1[0]*E1[0]+E1[2]*E1[2])*rechdir
      if Eabs == 0: return
      r1 = r + E1*d/2
      E2 = Edipole(r1,t)
      r2 = r + E2*d/2
      E3 = Edipole(r2,t)
      r3 = r + E3*d
      E4 = Edipole(r3,t)
      rn = r+(E1/6+E2/3+E3/3+E4/6)*d
      return rn
    
    # draw upper part of lobe
    x1 = []
    z1 = []
    x1.append(rstart[0])
    z1.append(0)
    rechdir = 1
    rpt = 0
    r = rstart
    rn = rnext(r,t,rechdir)
    rpt = 0
    nullct = 0
    # limit drawing not closer than 0.05 and no larger steps than 0.2
    # count zero passing, after one passing, stop
    while rpt<5000 and np.linalg.norm(r)>mindistance and np.linalg.norm(r-rn)<maxstep:        
       # next point 
       rn = rnext(r,t,rechdir)
       if rn[2]>0 and r[2]<0: nullct += 1 # from bottom to top
       if rn[2]<0 and r[2]>0: nullct += 1 # from top to bottom
       if nullct == 1:
           x1.append(rn[0])
           z1.append(rn[2])           
           break  
       rpt += 1
       r = rn
       x1.append(rn[0])
       z1.append(rn[2])
    
    # draw lower part of lobe   
    # count zero passing, after one passing, stop
    x2 = []
    z2 = []
    x2.append(rstart[0])
    z2.append(0)
    rechdir = -1
    rpt = 0
    r = rstart
    rn = rnext(r,t,rechdir)
    rpt = 0
    nullct = 0
    # limit drawing not closer than 0.05 and no larger steps than 0.2
    while rpt<5000 and np.linalg.norm(r)>mindistance and np.linalg.norm(r-rn)<maxstep:        
       # next point 
       rn = rnext(r,t,rechdir)
       if rn[2]>0 and r[2]<0: nullct += 1 # from bottom to top
       if rn[2]<0 and r[2]>0: nullct += 1 # from top to bottom
       if nullct == 1:
           x2.append(rn[0])
           z2.append(rn[2])           
           break  
       rpt += 1
       r = rn
       x2.append(rn[0])
       z2.append(rn[2])
                  
    return np.array([x1,z1,x2,z2])


# --------------------------------------------------------------------------
# setup plot
# -------------------------------------------------------------------------
ymax = xmax
fig, ax = plt.subplots(1,1,figsize=(8,4.5))
ax.set(xlabel="x (m)",ylabel="y (m)")
ax.set(xlim=(-xmax,xmax),ylim=(-ymax,ymax))

poynting = np.zeros((nbpoynting,nbpoynting))

# field E**2
im = plt.imshow(poynting,extent=(-xmax,xmax,-ymax,ymax),cmap='RdBu',label='abs(E)',alpha=0.5)
fig.colorbar(im,ax=ax,label='abs(E) (a.u.)')

# info box
infobox = ''
infobox += 'Wavelength: '+ "{:.1f}".format(lamb)+' (m)\n'
infobox += 'Frequency: '+ "{:.1e}".format(f)+' (1/s)'
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.98,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)   

scalearrows = 20
lobelw = 0.5

dtopx = [0,0]
dtopy = [0,0]
dbottomx = [0,0]
dbottomy = [0,0]
dipoletop = ax.plot(dtopx,dtopy,color='r',lw=1,alpha=0.8, label ='')[0]           
dipolebottom = ax.plot(dbottomx,dbottomy,color='b',lw=1,alpha=0.8, label ='')[0]           
dipoletopmarker = ax.plot(0,0.5,color='r',markersize=2,lw=0,marker='o',alpha=0.8, label ='')[0]           
dipolebottommarker = ax.plot(0,-0.5,color='b',markersize=2,lw=0,marker='o',alpha=0.8, label ='')[0]           

# --------------------------------------------------------------------------
# Main loop
# -------------------------------------------------------------------------

lobeleft = []
xstartleft = []
quiverleft = []

loberight = []
xstartright = []
quiverright = []

nblobes = 0
t = 0
time = 0
print('New lobe every: ',"{:.2e}".format(dtnewlobe)+' ns')

def animate(k):
    
   global t,time,lobeleft,loberight,xstartleft,xstartright,nblobes
  
   t += dt 
   time += dt
   # propagate existing lobes
   for k in range(0,nblobes,2):
       xstartright[k] += c*dt
       xstart = xstartright[k] 
       rplot = np.array([xstart,0,0])
       x1,z1,x2,z2 = Efieldline2(rplot,t)
       loberight[k].set_xdata(x1)
       loberight[k].set_ydata(z1)
       loberight[k+1].set_xdata(x2)
       loberight[k+1].set_ydata(z2)
       quiverright[k].set_offsets([xstart,0])
       quiverright[k+1].set_offsets([xstart,0])

       xstartleft[k] -= c*dt
       xstart = xstartleft[k] 
       rplot = np.array([xstart,0,0])
       x1,z1,x2,z2 = Efieldline2(rplot,t)
       lobeleft[k].set_xdata(x1)
       lobeleft[k].set_ydata(z1)
       lobeleft[k+1].set_xdata(x2)
       lobeleft[k+1].set_ydata(z2)
       quiverleft[k].set_offsets([xstart,0])
       quiverleft[k+1].set_offsets([xstart,0])
    
   # draw central dipoole animation    
   px,py,pz = dipole(t)
   ptop = 0.5/p0*pz
   pbottom = -ptop
   dipoletopmarker.set_xdata([0])
   dipoletopmarker.set_ydata([ptop])
   dipolebottommarker.set_xdata([0])
   dipolebottommarker.set_ydata([pbottom])
       
   dtopx = [0,0]
   dtopy = [0,ptop]
   dbottomx = [0,0]
   dbottomy = [0,pbottom]
   dipoletop.set_xdata(dtopx)
   dipoletop.set_ydata(dtopy)
   dipolebottom.set_xdata(dbottomx)
   dipolebottom.set_ydata(dbottomy)

   # Update E**2 plot
   for i in range(nbpoynting):
      for j in range(nbpoynting):
          xp = -xmax+i/nbpoynting*2*xmax 
          zp = -ymax+j/nbpoynting*2*xmax
          rp = np.array([xp,0,zp])
          Ep = Edipole(rp,t)
          #  Hp = Hdipole(rp,t)
          #  poynting[i,j] = np.linalg.norm(np.cross(Ep,Hp))/poyntingscale
          if Ep[2]>0:
             poynting[i,j] = np.linalg.norm(Ep)/poyntingscale
          else:  
             poynting[i,j] = -np.linalg.norm(Ep)/poyntingscale
   im.set_array(np.transpose(poynting))       
       
   # for each dtnewlobe create new  at origin
   if time>=dtnewlobe: 
       xstart = 0
       rplot = np.array([xstart,0,0])
       x1,z1,x2,z2 = Efieldline2(rplot,t)
       of = orientationfieldline(t)
       if of == +1: # -dtnewlobe
           xstartright.append(xstart)
           xstartright.append(xstart) 
           loberight.append(ax.plot(x1,z1,color='b',lw=lobelw,alpha=0.8, label ='')[0])             
           loberight.append(ax.plot(x2,z2,color='b',lw=lobelw,alpha=0.8, label ='')[0]) 
           quiverright.append(ax.quiver(0, 0, 0, 1, pivot = 'middle', scale = scalearrows, color="b", alpha=0.5))
           quiverright.append(ax.quiver(0, 0, 0, 1, pivot = 'middle', scale = scalearrows, color="b", alpha=0.5))           
        
           xstartleft.append(xstart)
           xstartleft.append(xstart) 
           lobeleft.append(ax.plot(x1,z1,color='b',lw=lobelw,alpha=0.8, label ='')[0])             
           lobeleft.append(ax.plot(x2,z2,color='b',lw=lobelw,alpha=0.8, label ='')[0])             
           quiverleft.append(ax.quiver(0, 0, 0, 1, pivot = 'middle', scale = scalearrows, color="b", alpha=0.5))
           quiverleft.append(ax.quiver(0, 0, 0, 1, pivot = 'middle', scale = scalearrows, color="b", alpha=0.5))           
           nblobes += 2
       elif of == -1:
           xstartright.append(xstart)
           xstartright.append(xstart) 
           loberight.append(ax.plot(x1,z1,color='r',lw=lobelw,alpha=0.8, label ='')[0])                       
           loberight.append(ax.plot(x2,z2,color='r',lw=lobelw,alpha=0.8, label ='')[0])                       
           quiverright.append(ax.quiver(0, 0, 0, -1, pivot = 'middle', scale = scalearrows, color="r", alpha=0.5))
           quiverright.append(ax.quiver(0, 0, 0, -1, pivot = 'middle', scale = scalearrows, color="r", alpha=0.5))           
           
           xstartleft.append(xstart)
           xstartleft.append(xstart) 
           lobeleft.append(ax.plot(x1,z1,color='r',lw=lobelw,alpha=0.8, label ='')[0])             
           lobeleft.append(ax.plot(x2,z2,color='r',lw=lobelw,alpha=0.8, label ='')[0])             
           quiverleft.append(ax.quiver(0, 0, 0, -1, pivot = 'middle', scale = scalearrows, color="r", alpha=0.5))
           quiverleft.append(ax.quiver(0, 0, 0, -1, pivot = 'middle', scale = scalearrows, color="r", alpha=0.5))           
           nblobes += 2

       
       time = 0
       print('Lobe: ',nblobes,' at ',"{:.2e}".format(t)+' ns of',' at ',"{:.2e}".format(tend)+' ns')
  
   fig.suptitle('time: '+"{:.3f}".format(t/1e-9)+' ns')
     
    
anim = animation.FuncAnimation(fig,animate,interval=1,frames=int(tend/dt))
anim.save(animationfile,fps=25,dpi=300)

