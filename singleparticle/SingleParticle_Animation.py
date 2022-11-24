#
# Single particle movement
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

# model = "Gyration"
# model = "ExB"
model = "Dipole"
# model = "Gradientx"


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
if model == "Gyration":
    Bz0 = 0.05 # B field in z direction
    Ex0 = 0 # E field in x direction
    nu = 2e8  #collision frequency for e scattering
    animationname = 'ParticleGyration.gif'
    Te = 3 # Elektronentemperatur in eV
    v0e = np.sqrt(Te*np.abs(el)/me) # thermal velocity electron
    v0z = v0e/10 # initial z velocity
    x0z = -2e-4 # x initial z position
    elev = 30 # azimuthal viewing angle
    sam = 5e-4 # +- simulation box size in m
    scaleb=50 # scale factor B field arrows
    scalee=1 # scale factor E field arrows
elif model =="ExB":
    Bz0 = 0.05 # B field in z direction
    Ex0 = 5000 # E field in x direction
    nu = 2e8 # collision frequency for e scattering
    animationname = 'ParticleExB.gif'  
    Te = 3 # Elektronentemperatur in eV
    v0e = np.sqrt(Te*np.abs(el)/me) # thermal velocity electron
    v0z = v0e/10 # initial z velocity   
    x0z = -2e-4 # x initial z position
    elev = 30 # azimuthal viewing angle
    sam = 5e-4 # +- simulation box size in m
    scaleb=30 # scale factor B field arrows
    scalee=3e-9 # scale factor E field arrows
elif model == "Dipole":
    Bz0 = 0.1 # B field in z direction
    Ex0 = 0 # E field in x direction
    nu = 0 # collision frequency for e scattering
    Te = 3 # Elektronentemperatur in eV
    animationname = 'ParticleBottle.mp4'    
    v0e = np.sqrt(Te*np.abs(el)/me) # thermal velocity electron
    v0z = v0e/50 # initial z velocity
    x0z = -2e-4 # x initial z position
    mu = 2.1e-5 * np.array([0.0, 0.0, 1.0]) # set magnetic moment to point in z direction
    elev = 5 # azimuthal viewing angle
    sam = 5e-4 # +- simulation box size in m
    z_disp = 5e-4 # displacement of the two magnetic dipoles with respect to zero (one at z = -z_disp, the other at +z_disp)
    scaleb=20 # scale factor B field arrows
    scalee=1 # scale factor E field arrows
elif model == "Gradientx":
    Bz0 = 0.1 # B field in z direction
    Ex0 = 0 # E field in x direction
    Te = 3 # Elektronentemperatur in eV
    nu = 0 # collision frequency for e scattering
    animationname = 'ParticleGradientx.gif'
    v0e = np.sqrt(Te*np.abs(el)/me) # thermal velocity electron
    v0z = 0 # initial z velocity
    x0z = -2e-4 # x initial z position
    elev = 30 # azimuthal viewing angle
    sam = 5e-4 # +- simulation box size in m
    scaleb=5 # scale factor B field arrows
    scalee=1 # scale factor E field arrows
elif model == "Magnetron":
    Bz0 = 0.1 # B field in z direction
    Ex0 = 0 # E field in x direction
    nu = 0 # collision frequency for e scattering
    Te = 3 # Elektronentemperatur in eV
    animationname = 'ParticleMagnetron.mp4'    
    v0e = np.sqrt(Te*np.abs(el)/me) # thermal velocity electron
    v0z = v0e/50 # initial z velocity
    x0z = -2e-4 # x initial z position
    mu = 2.1e-5 * np.array([0.0, 0.0, 1.0]) # set magnetic moment to point in z direction
    elev = 5 # azimuthal viewing angle
    sam = 30e-3 # +- simulation box size in m
    z_disp = 5e-4 # displacement of the two magnetic dipoles with respect to zero (one at z = -z_disp, the other at +z_disp)
    scaleb=20 # scale factor B field arrows
    scalee=1 # scale factor E field arrows    
    

# --------------------------------------------------------------------------
# Define fields
# --------------------------------------------------------------------------
def Bfield(x):
    global model, z_displ
    if model == "Gyration":        
        Bx = 0
        By = 0
        Bz = Bz0 
    elif model == "ExB":
        Bx = 0
        By = 0
        Bz = Bz0        
    elif model =="Gradientx":
        Bx = 0
        By = 0
        Bz = (1+x[0]*2e3)*Bz0
    elif model == "Dipole":
        # point dipole A
        pos_A = np.array([0.0, 0.0, z_disp]) # first dipole
        r_A = x - pos_A                      
        rmag_A = np.sqrt(sum(r_A**2))
        B1_A = 3.0*r_A*np.dot(mu,r_A) / (rmag_A**5)   # calculate the first term to the magnetic field
        B2_A = -1.0 * mu / (rmag_A**3)                # calculate the second term
    
        # point dipole B
        pos_B = np.array([0.0, 0.0, -z_disp])  # set the position of the first dipole
        r_B = x - pos_B        # find the difference between this position and the observation position
        rmag_B = np.sqrt(sum(r_B**2))
        B1_B = 3.0*r_B*np.dot(mu,r_B) / (rmag_B**5) # calculate the first term to the magnetic field
        B2_B = -1.0 * mu / (rmag_B**3)              # calculate the second term
    
        Bx =  ((mu_0/(4.0*np.pi)) * (B1_A + B2_A + B1_B + B2_B))[0]  
        By =  ((mu_0/(4.0*np.pi)) * (B1_A + B2_A + B1_B + B2_B))[1]  
        Bz =  ((mu_0/(4.0*np.pi)) * (B1_A + B2_A + B1_B + B2_B))[2]       
        
    return np.array([Bx,By,Bz]) # return field due to magnetic bottle


def Efield(x):
    Ex = Ex0
    Ey = 0
    Ez = 0
    return np.array([Ex,Ey,Ez])




# Initialize Electron Position, Velocity, Acceleration
v = np.array([v0e,0,v0z])
x = np.array([0,0,x0z])
a = np.array([0,0,0])

# reserve memory for trajectorys
xt = np.array(0)
yt = np.array(0)
zt = np.array(-2e-4)

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
            if model == 'dipole' or model == 'magnetron': # choose simple arrow lengths for model dipole
                lengthb[i,j,m] = 0.2/scaleb 
            axp.quiver(xgrid[i,j,m], ygrid[i,j,m], zgrid[i,j,m], 
                       ubfeld[i,j,m], vbfeld[i,j,m], wbfeld[i,j,m], 
                       pivot = 'tail', color="b", length=lengthb[i,j,m]*scaleb, alpha=0.5)


uefeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
vefeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
wefeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
lengthe = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
for i in range(len(xgrid[:,0,0])):
    for j in range(len(xgrid[0,:,0])):
        for m in range(len(xgrid[0,0,:])):
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
if model == "Gyration":
    infobox += 'B: ' + "{:.2f}".format(Bz0) + ' (T)\n' 
elif model =="ExB":
    infobox += 'B: ' + "{:.2f}".format(Bz0) + ' (T)\n' 
    infobox += 'E: ' + "{:.2f}".format(Ex0) + ' (V/m)\n' 
elif model =="Dipole":
    infobox += 'B: dipole\n'     
elif model =="Gradientx":
    infobox += 'B: gradient in x\n' 
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
   
        
    t += dt
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
  
         
       
# ----------------------------------------
# Create the complete time sequence
# teh total time span is tend/dtplot
# returns the rate in
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=180)

