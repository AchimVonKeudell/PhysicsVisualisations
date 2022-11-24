#
# Single particle movement in a 2'' magnetron
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Implementation of the Bfield by
# Julian Held
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

#model = 'GradientBDrift'
#model = 'ExBDrift'
model = 'OpenFieldLine'

if model == 'GradientBDrift':
    Ex0 = 0 # E field in x direction
    Ez0 = 0
    nu = 0 # collision frequency for e scattering
    x0z = 1e-3 # x initial z position
    x0x = 15e-3 # x initial x position
    animationname = 'ParticleMagnetron_GradB.mp4'    
if model == 'GradientBDriftColl':
    Ex0 = 0 # E field in x direction
    Ez0 = 0
    nu = 1e7 # collision frequency for e scattering
    x0z = 1e-3 # x initial z position
    x0x = 15e-3 # x initial x position
    animationname = 'ParticleMagnetron_GradB_coll.mp4'    
if model == 'ExBDrift':
    Ex0 = 0 # E field in x direction
    Ez0 = -100
    nu = 0 # collision frequency for e scattering    
    x0z = 1e-3 # x initial z position
    x0x = 15e-3 # x initial x position
    animationname = 'ParticleMagnetron_ExB.mp4'    
if model == 'OpenFieldLine':
    Ex0 = 0 # E field in x direction
    Ez0 = 0
    nu = 0 # collision frequency for e scattering
    x0z = 1e-3 # x initial z position
    # x0x = 16e-3 # x stays inside
    # x0x = 17e-3 # almost deconfinement
    # x0x = 18e-3 # quick deconfinement
    # x0x = 19e-3 # almost confined
    x0x = 20e-3 # almost confined
    animationname = 'ParticleMagnetron_OpenField_x20mm.mp4'    

    
# --------------------------------------------------------------------------
# Timing Simulation
# --------------------------------------------------------------------------
t = 0
dt = 1e-11 # time step
dtplot = 2e-10
time = 0
tend = 2e-7
steps = int(tend/dtplot)
trotation = 4e-7 # time for totation of the image around z-axis


# --------------------------------------------------------------------------
# Physical parameter Simulation
# --------------------------------------------------------------------------

Te = 3 # Elektronentemperatur in eV
v0e = np.sqrt(Te*np.abs(el)/me) # thermal velocity electron
v0z = v0e*10 # initial z velocity
elev = 25 # azimuthal viewing angle
sam = 30e-3 # +- simulation box size in m
scaleb=500 # scale factor B field arrows
scalee=0.0003 # scale factor E field arrows    
    

# --------------------------------------------------------------------------
# Define fields
# --------------------------------------------------------------------------
def Bfield(x):
    """ Returns Br, Bz in T for r, z position (in m).
    Reconstructed for the Thin Films consulting 2" UHV magnetrons.
    Reconstruction performed by Dennis Krueger (TET).
    z = 0 is the target surface, if traget_thickness is set correctly.
    z = 0 is the magnetron surface if target_thickness=0

    If used, cite: Krueger et al, Phys. Plasmas 25, 061207 (2018)"""
    target_thickness=3e-3
    R = 20.57870655967374e-3 #radius of the ring magnet (m)
    dC = 16.68888450585952e-3 #location of center magnet below the target (m)
    dR = 8.240298921866946e-3 #location of the ring magnet below the target (m)
    mC = -6.07855 #magnetic dipole moment of the center magnet (Am^2)
    mR = 10.2386 #magnetic dipole moment of the ring magnet (Am^2)

    r = np.sqrt(x[0]**2+x[1]**2)
    z = x[2]
    z = z + target_thickness

    Br = (3. * const.mu_0 * mC/(4.*np.pi) * (r*(z + dC)) / (((z+dC)**2 + r**2)**(5./2.))
    + 3. * const.mu_0 * mR / (4.*np.pi)  * (r*(z + dR))
    / (((z+dR)**2 + r**2 +R**2)**(5./2.))
    * hyp2f1(3./4., 5./4., 1., (4. * r**2 * R**2)
    / (((z+dR)**2 + r**2 + R**2)**2))
    - 15. * const.mu_0 * mR / (8.*np.pi) *  (r*R**2*(z+dR) *( (z+dR)**2-r**2+R**2)
    /(((z+dR)**2 + r**2 + R**2))**(9./2.))
    * hyp2f1(7./4., 9./4., 2., (4. * r**2 * R**2)
    / (((z+dR)**2+r**2+R**2)**2)))


    Bz = (const.mu_0 * mC/(4.*np.pi) * (2.*(z+dC)**2 - r**2)/((z+dC)**2 + r**2)**(5./2.)
       + const.mu_0 * mR/(4.*np.pi) * (2.*(z+dR)**2 - r**2 - R**2)
       / ((z+dR)**2 + r**2 + R**2)**(5./2.)
       * hyp2f1(3./4., 5./4., 1., (4. * r**2 * R**2)
       / ( ((z+dR)**2 + r**2 + R**2)**2) )
       + 15. * const.mu_0 * mR / (4.*np.pi) * (r**2*R**2*(z+dR)**2)
       / (((z+dR)**2 + r**2 + R**2)**(9./2.))
       * hyp2f1(7./4., 9./4., 2., (4. * r**2 * R**2)/((z+dR)**2+r**2+R**2)**2))

    if x[0] == 0:
       if x[1]>=0:
          phi = np.pi/2
       else:
          phi = 3/2*np.pi
    else:
       if x[0] >=0 and x[1]>=0: 
         phi = np.arctan(abs(x[1]/x[0]))
       elif x[0] <=0 and x[1]>=0: 
         phi = np.pi-np.arctan(abs(x[1]/x[0]))           
       elif x[0] <=0 and x[1]<=0: 
         phi = np.pi+np.arctan(abs(x[1]/x[0]))           
       elif x[0]>=0 and x[1]<=0: 
         phi = 2*np.pi-np.arctan(abs(x[1]/x[0]))           
         
    
    Bx = Br*np.cos(phi)
    By = Br*np.sin(phi)        
        
    return np.array([Bx,By,Bz]) # return field due to magnetic bottle


def Efield(x):
    Ex = Ex0
    Ey = 0
    Ez = Ez0
    return np.array([Ex,Ey,Ez])




# Initialize Electron Position, Velocity, Acceleration
v = np.array([v0e,0,v0z])
x = np.array([x0x,0,x0z])
a = np.array([0,0,0])

# reserve memory for trajectorys
xt = np.array(x0x)
yt = np.array(0)
zt = np.array(x0z)

# --------------------------------------------------------------------------
# Setup Plot
# --------------------------------------------------------------------------
sa = sam/1e-3 # plot axis scale in mm
hvsw = 1
gxy = 10
gz = 10

fig = plt.figure(figsize=(8,4.5))
axp = plt.axes(projection='3d')
axp.set(xlabel="x (mm)",ylabel="y (mm)",zlabel="z (mm)")
axp.set(xlim=(-sa,sa),ylim=(-sa,sa),zlim=(0,2*hvsw*sa))

# trajectory line
line1 = axp.plot3D(xt,yt,zt,color='r',lw=1)[0]

# marker at the end
xp = np.array(1)
yp = np.array(1)
zp = np.array(1)
line1p = axp.plot3D(xp,yp,zp,color='r',marker='o',markersize=2,lw=0)[0]


# define field arrows in the plot in mm
# quiver can only generate arrows with different length by
# using individual axp.quiver commands
xgrid, ygrid, zgrid = np.meshgrid(np.arange(-sa,sa+sa/(gxy),2*sa/(gxy)),
                                  np.arange(-sa,sa+sa/(gxy),2*sa/(gxy)),
                                  np.arange(0,hvsw*2*sa+hvsw*sa/gz,hvsw*2*sa/gz))

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
            axp.quiver(xgrid[i,j,m], ygrid[i,j,m], zgrid[i,j,m], 
                       ubfeld[i,j,m], vbfeld[i,j,m], wbfeld[i,j,m], 
                       pivot = 'tail', color="b", length=lengthb[i,j,m]*scaleb, alpha=0.5)


uefeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
vefeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
wefeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
lengthe = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
for i in range(len(xgrid[:,0,0])):
    for j in range(len(xgrid[0,:,0])):
            m = 3
            # real points in m, xgrid points in mm
            xb = xgrid[i,j,m]*1e-3
            yb = ygrid[i,j,m]*1e-3
            zb = zgrid[i,j,m]*1e-3
            uefeld[i,j,m], vefeld[i,j,m], wefeld[i,j,m] = Efield([xb,yb,zb])
            lengthe[i,j,m] = np.sqrt(Efield([xb,yb,zb])[0]**2+Efield([xb,yb,zb])[1]**2+Efield([xb,yb,zb])[2]**2)
            axp.quiver(xgrid[i,j,m], ygrid[i,j,m], zgrid[i,j,m], 
                       uefeld[i,j,m], vefeld[i,j,m], wefeld[i,j,m], 
                       pivot= 'tail', color="g", length=lengthe[i,j,m]*scalee, alpha = 0.5)
            
# plotting the arrows at 0,0,0 twice to generate legend            
axp.quiver(xgrid[0,0,0], ygrid[0,0,0], zgrid[0,0,0], 
                       ubfeld[0,0,0], vbfeld[0,0,0], wbfeld[0,0,0], 
                       pivot = 'tail', length=lengthb[0,0,0]*scaleb,color="b",label="B field", alpha = 0.5)
axp.quiver(xgrid[i,j,m], ygrid[0,0,0], zgrid[0,0,0], 
                       uefeld[0,0,0], vefeld[0,0,0], wefeld[0,0,0], 
                       pivot = 'tail', length=lengthe[0,0,0]*scalee,color="g",label="E field", alpha = 0.5)
axp.legend(fontsize=6)

# info box
infobox = ''
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
  
  time = 0
  while time < dtplot:    
    t += dt
    time += dt
    # ------------------------------------------------
    # Boris-Bunemann Pusher
    # ------------------------------------------------
        
    ## Adding 1st half 
    v += el/me * Efield(x) * 0.5*dt 

    ## Rotation according to B
    tw = el/me * dt * Bfield(x) * 0.5 
    v1 = v + np.cross(v, tw)
    v += np.cross(v1, 2 /(1 + np.abs(tw)**2) * tw)

    ## Adding 2nd half
    v += el/me * Efield(x) * dt * 0.5

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
    
  # rotate the view
  axp.view_init(elev,45+t/trotation*360)
    
  # Update the title of the plot
  fig.suptitle('time: ' + "{:.2e}".format(t)+' s')
  print('time: ' + "{:.2e}".format(t)+' s of '+ "{:.2e}".format(tend)+' s')
           
       
# ----------------------------------------
# Create the complete time sequence
# the total time span is tend/dtplot
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=180)

