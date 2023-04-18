#
# Single particle movement in a cyclotron
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as const
from scipy.special import hyp2f1

el = -1.6e-19 # el charge
me = 9e-31 # el mass
mp = 1.6e-27 # proton mass
mu_0 = np.pi * 4.0e-7 # permeability of free space [kg*m*s^-2*A^-2]

# --------------------------------------------------------------------------
# Choose Model 
# --------------------------------------------------------------------------

model = 'Cyclotron'


# --------------------------------------------------------------------------
# Timing Simulation
# --------------------------------------------------------------------------
t = 0
dt = 1e-11 # time step
steps = 1000 # simulation steps
trotation = 1e-8 # time for totation of the image around z-axis


# --------------------------------------------------------------------------
# Physical parameter Simulation
# --------------------------------------------------------------------------
if model == "Cyclotron":
    Bz0 = 0.05 # B field in z direction
    Ex0 = 5000 # E field in x direction
    nu = 0  #collision frequency for e scattering
    animationname = 'ParticleCyclotron.mp4'
    Te = 3 # Elektronentemperatur in eV
    v0e = np.sqrt(Te*np.abs(el)/me) # thermal velocity electron
    v0z = 0 # initial z velocity
    x0y = -7e-5
    x0z = -1e-6 # x initial z position
    elev = 30 # azimuthal viewing angle
    sam = 5e-4 # +- simulation box size in m
    scaleb=50 # scale factor B field arrows
    scalee=3e-9 # scale factor E field arrows
    omegaE = el*Bz0/me # cyclotron frequency
    rcyclotron = 0.0004 # radius in m
    
    

# --------------------------------------------------------------------------
# Define fields
# --------------------------------------------------------------------------
def Bfield(x):
    global model, z_displ

    Bx = 0
    By = 0
    if x[0]**2+x[1]**2<rcyclotron**2:
      Bz = Bz0
    else:
      Bz = 0
        
    return np.array([Bx,By,Bz]) # return field due to magnetic bottle


def Efield(x,t):
      Ex = 0
      Ey = 0
      Ez = 0
      if ((x[0]>-0.001) and (x[0]<0.001)):
        Ex = Ex0*np.cos(omegaE*t)    
    
      return np.array([Ex,Ey,Ez])



# Initialize Electron Position, Velocity, Acceleration
v = np.array([v0e,0,v0z])
x = np.array([0,x0y,x0z])
a = np.array([0,0,0])

# reserve memory for trajectorys
xt = np.array(0)
yt = np.array(x0y)
zt = np.array(x0z)

# --------------------------------------------------------------------------
# Setup Plot
# --------------------------------------------------------------------------
sa = sam/1e-3 # plot axis scale in mm
hvsw = 1
gxy = 5
gz = 5

fig = plt.figure(figsize=(8,4.5))
axp = plt.axes(projection='3d')
axp.set(xlabel="x (mm)",ylabel="y (mm)",zlabel="z (mm)")
axp.set(xlim=(-sa,sa),ylim=(-sa,sa),zlim=(-hvsw*sa,hvsw*sa))

# trajectory line
line1 = axp.plot3D(xt,yt,zt,color='r',lw=1)[0]

# marker at the end
xp = np.array(1)
yp = np.array(1)
zp = np.array(1)
line1p = axp.plot3D(xp,yp,zp,color='r',marker='o',markersize=3,lw=0)[0]


# define field arrows in the plot in mm
# quiver can only generate arrows with different length by
# using individual axp.quiver commands
xgrid, ygrid, zgrid = np.meshgrid(np.arange(-sa,sa+2*sa/(gxy),2*sa/(gxy)),
                                  np.arange(-sa,sa+2*sa/(gxy),2*sa/(gxy)),
                                  np.arange(-sa*hvsw,hvsw*sa+2*hvsw*sa/gz,2*hvsw*sa/gz))

ubfeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
vbfeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
wbfeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
lengthb = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
for i in range(len(xgrid[:,0,0])):
    for j in range(len(xgrid[0,:,0])):
        for m in range(len(xgrid[0,0,:])):
            # real points in m, xgrid points in mm
            xb = xgrid[i,j,m]*1e-3
            yb = ygrid[i,j,m]*1e-3
            zb = zgrid[i,j,m]*1e-3
            ubfeld[i,j,m], vbfeld[i,j,m], wbfeld[i,j,m] = Bfield([xb,yb,zb])
            lengthb[i,j,m] = np.sqrt( Bfield([xb,yb,zb])[0]**2 + 
                                      Bfield([xb,yb,zb])[1]**2 + Bfield([xb,yb,zb])[2]**2)
            
            if xb**2+yb**2<rcyclotron**2:
              axp.quiver(xgrid[i,j,m], ygrid[i,j,m], zgrid[i,j,m], 
                         ubfeld[i,j,m], vbfeld[i,j,m], wbfeld[i,j,m], 
                         pivot = 'tail', color="b", length=lengthb[i,j,m]*scaleb, alpha=0.5)


# plotting e-field
quiverefeld = []
yefeld = np.arange(-sa,sa+2*sa/(gxy),2*sa/(gxy))
nbquiver = len(yefeld)
xefeld = np.zeros(nbquiver)
zefeld = np.zeros(nbquiver)               
for i in range(nbquiver):
        quiverefeld.append(axp.quiver(xefeld[i], yefeld[i], zefeld[i], 1, 0, 0, pivot = 'middle', length = 0.1, color="g", alpha=0.5))
  
quiverefeld1 = axp.quiver(xefeld[0], yefeld[0], zefeld[0], 0, 1, 0, pivot = 'middle', length = 0.0, color="g", alpha=0.5,label='E field')        

# plotting d shaped magnet
angle = np.linspace(np.pi/2+np.pi*0.05,np.pi/2+np.pi*0.95,180)
xdmagnet = rcyclotron*np.cos(angle)/1e-3      
ydmagnet = rcyclotron*np.sin(angle)/1e-3
zdmagnet = 1e-4*np.ones(180)/1e-3
xdmagnet = np.append(xdmagnet,xdmagnet[0])
ydmagnet = np.append(ydmagnet,ydmagnet[0])
zdmagnet = np.append(zdmagnet,zdmagnet[0])
axp.plot3D(xdmagnet,ydmagnet,zdmagnet,color='orange',lw=3)      

angle = np.linspace(np.pi/2+np.pi*0.05,np.pi/2+np.pi*0.95,180)
xdmagnet = rcyclotron*np.cos(angle)/1e-3      
ydmagnet = rcyclotron*np.sin(angle)/1e-3
zdmagnet = -1e-4*np.ones(180)/1e-3
xdmagnet = np.append(xdmagnet,xdmagnet[0])
ydmagnet = np.append(ydmagnet,ydmagnet[0])
zdmagnet = np.append(zdmagnet,zdmagnet[0])
axp.plot3D(xdmagnet,ydmagnet,zdmagnet,color='orange',lw=3)      

angle = np.linspace(np.pi/2+np.pi*1.05,np.pi/2+np.pi*1.95,180)
xdmagnet = rcyclotron*np.cos(angle)/1e-3      
ydmagnet = rcyclotron*np.sin(angle)/1e-3
zdmagnet = 1e-4*np.ones(180)/1e-3
xdmagnet = np.append(xdmagnet,xdmagnet[0])
ydmagnet = np.append(ydmagnet,ydmagnet[0])
zdmagnet = np.append(zdmagnet,zdmagnet[0])
axp.plot3D(xdmagnet,ydmagnet,zdmagnet,color='orange',lw=3)      
            
angle = np.linspace(np.pi/2+np.pi*1.05,np.pi/2+np.pi*1.95,180)
xdmagnet = rcyclotron*np.cos(angle)/1e-3      
ydmagnet = rcyclotron*np.sin(angle)/1e-3
zdmagnet = -1e-4*np.ones(180)/1e-3
xdmagnet = np.append(xdmagnet,xdmagnet[0])
ydmagnet = np.append(ydmagnet,ydmagnet[0])
zdmagnet = np.append(zdmagnet,zdmagnet[0])
axp.plot3D(xdmagnet,ydmagnet,zdmagnet,color='orange',lw=3)      

# plotting the arrows at 0,0,0 twice to generate legend            
axp.quiver(xgrid[0,0,0], ygrid[0,0,0], zgrid[0,0,0], 
                       ubfeld[0,0,0], vbfeld[0,0,0], wbfeld[0,0,0], 
                       pivot = 'tail', length=lengthb[0,0,0]*scaleb,color="b",label="B field", alpha = 0.5)
axp.legend(fontsize=6)

# info box
infobox = ''
infobox += 'B: ' + "{:.2f}".format(Bz0) + ' (T)\n' 
infobox += 'w_ce: ' + "{:.2e}".format(el*Bz0/me) + ' (1/s)\n' 
infobox += 'nu: ' + "{:.0e}".format(nu) + ' (1/s)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
axp.text2D(0,1,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=axp.transAxes)



# -----------------------------------------------------------------------------
# Collision calculation electrons
#  from eduPIC code example by Z. Donko et al.
# -----------------------------------------------------------------------------
def collision(v):
        vb = np.sqrt(v[0]**2+v[1]**2+v[2]**2)
        # calulate Euler angles
        if v[0] == 0:
            theta = 0.5*np.pi
        else:
            theta = np.arctan(np.sqrt(v[1]**2+v[2]**2)/v[0])
        if v[1] == 0:
            if v[2] > 0:
                phi = 0.5*np.pi
            else:
                phi = -0.5*np.pi
        else:        
            phi = np.arctan(v[2]/v[1])

        st = np.sin(theta)
        ct = np.cos(theta)
        sp = np.sin(phi)
        cp = np.cos(phi)

        # select scattering angles 
        xi = np.arccos(1-2*np.random.rand())
        eta = 2*np.pi*np.random.rand()
        
        sc = np.sin(xi)
        cc = np.cos(xi)
        se = np.sin(eta)
        ce = np.cos(eta)
        
        # rotate velocity vector
        vx = vb * (ct * cc - st * sc * ce);
        vy = vb * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
        vz = vb * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
   
        return np.array([vx,vy,vz]) 
   
   
# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
    global t,v,x,a
    global xt,yt,zt
    global quiverefeld
   
    print(k,' of ',steps)    
    
    t += dt
    # ------------------------------------------------
    # Boris-Bunemann Pusher
    # ------------------------------------------------
        
    ## Adding 1st half 
    v += el/me * Efield(x,t) * 0.5*dt 

    ## Rotation according to B
    tw = el/me * dt * Bfield(x) * 0.5 
    v1 = v + np.cross(v, tw)
    v += np.cross(v1, 2 /(1 + np.abs(tw)**2) * tw)

    ## Adding 2nd half
    v += el/me * Efield(x,t) * dt * 0.5   

    # Updating position
    x += dt * v      
    
    # check after microsteps*dt whether a collsion should take place
    if np.random.rand() < 1-np.exp(-nu*dt):
        v = collision(v)
    
    # Update the new positions in vector
    xt = np.append(xt,x[0])
    yt = np.append(yt,x[1])
    zt = np.append(zt,x[2])
    
    # real points in m, xgrid points in mm
    line1.set_xdata(xt/1e-3)
    line1.set_ydata(yt/1e-3) 
    line1.set_3d_properties(zt/1e-3) 
    xp = np.array([x[0]])
    yp = np.array([x[1]])
    zp = np.array([x[2]])
    line1p.set_xdata(xp/1e-3)
    line1p.set_ydata(yp/1e-3) 
    line1p.set_3d_properties(zp/1e-3) 

       #print(nbquiver) 
       #print('len1',len(quiverefeld))
       #for u in range(nbquiver):
       #   print(u) 
       #quiverefeld.clear
    if Efield(x,t)[0]>=0: 
        directionefield = 1
    else: 
        directionefield = -1
    elength = abs(Efield(x,t)[0])/Ex0*0.1    
    
    
    for i in range(nbquiver):
        quiverefeld[i].remove()
    #quiverefeld.clear
    quiverefeld = []   
    for i in range(nbquiver):
        quiverefeld.append(axp.quiver(xefeld[i], yefeld[i], zefeld[i], directionefield, 0, 0, pivot = 'middle', length = elength, color="g", alpha=0.5))
        
    # rotate the view
    axp.view_init(elev+t/trotation*(90-30),45+t/trotation*360)
    
    # Update the title of the plot
    fig.suptitle('time: ' + "{:.2e}".format(t)+' s')
  
         
       
# ----------------------------------------
# Create the complete time sequence
# teh total time span is tend/dtplot
# returns the rate in
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)

