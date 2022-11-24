#
# Single particle movement in a quadrupole
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as const
from scipy.special import hyp2f1


# --------------------------------------------------------------------------
# Choose Model 
# --------------------------------------------------------------------------
#model = 'massesclose'
model = 'massesfar'

el = 1.6e-19 # el charge
r0 = 0.0027

if model == 'massesclose':
    mion = 28
    mion2 = 27.95
    mp = mion*1.6e-27 # proton mass
    mp2 = mion2*1.6e-27 # proton mass
    dt = 1e-8 # time 
    f = 2e6
    tend = 50*1/f
    dtreplot = 2e-8
    animationname = 'qms_m28_m27.95.mp4'
    E0 = 3 # Ionen in eV
    v0e = np.sqrt(2*E0*np.abs(el)/mp) # thermal velocity electron
    v0z = v0e/10 # initial z velocity
    sam = 5e-3 # +- simulation box size in m
    z0 = -sam
    y0 = 3e-4
    x0 = 3e-4
    vmax = 8e4
    elev = 20 # azimuthal viewing angle
    scaleb=50 # scale factor B field arrows
    scalee=1 # scale factor E field arrows
    phi0 = 10
    U = 13
    V = 123.5
if model == 'massesfar':
    mion = 28
    mion2 = 27.50
    mp = mion*1.6e-27 # proton mass
    mp2 = mion2*1.6e-27 # proton mass
    dt = 1e-8 # time 
    f = 2e6
    tend = 50*1/f
    dtreplot = 2e-8
    animationname = 'qms_m28_m27.50.mp4'
    E0 = 3 # Ionen in eV
    v0e = np.sqrt(2*E0*np.abs(el)/mp) # thermal velocity electron
    v0z = v0e/10 # initial z velocity
    sam = 5e-3 # +- simulation box size in m
    z0 = -sam
    y0 = 3e-4
    x0 = 3e-4
    vmax = 8e4
    elev = 20 # azimuthal viewing angle
    scaleb=50 # scale factor B field arrows
    scalee=1 # scale factor E field arrows
    phi0 = 10
    U = 13
    V = 123.5



# --------------------------------------------------------------------------
# Timing Simulation
# --------------------------------------------------------------------------
t = 0
steps = int(tend/dtreplot) # simulation steps
trotation = tend # time for totation of the image around z-axis


# --------------------------------------------------------------------------
# gemetry parameter simulation
# --------------------------------------------------------------------------  
sa = sam/1e-3 # plot axis scale in mm
hvsw = 1
gxy = 5
gz = 5
    

# --------------------------------------------------------------------------
# Define fields
# --------------------------------------------------------------------------
def Efield(phiz,x):
    # phi = phi0 (x^2-y^2)/(2r0^2)
    Ex = phiz * 2 * x[0] /(2*r0**2)
    Ey = -phiz * 2 * x[1] /(2*r0**2)
    Ez = 0
    return np.array([Ex,Ey,Ez])

def phi(U,V,t):
    return U - V*np.cos(2*np.pi*f*t-np.pi*45/180)

# Initialize Electron Position, Velocity, Acceleration
x = np.array([x0,y0,z0])
v = np.array([0,0,v0z])
a = np.array([0,0,0])

# reserve memory for trajectorys
xt = np.array(x0)
yt = np.array(y0)
zt = np.array(z0)

# Initialize Electron Position, Velocity, Acceleration
x2 = np.array([x0,y0,z0])
v2 = np.array([0,0,v0z])
a2 = np.array([0,0,0])

# reserve memory for trajectorys
xt2 = np.array(x0)
yt2 = np.array(y0)
zt2 = np.array(z0)

# --------------------------------------------------------------------------
# Setup Plot
# --------------------------------------------------------------------------
fig = plt.figure(figsize=(8,4.5))
axp = fig.add_subplot(1,2,1,projection='3d')
axp.set(xlabel="x (mm)",ylabel="y (mm)",zlabel="z (mm)")
axp.set(xlim=(-sa,sa),ylim=(-sa,sa),zlim=(-hvsw*sa,hvsw*sa))

# trajectory line
line1 = axp.plot3D(xt,yt,zt,color='b',lw=1,label='ion '+"{:.2f}".format(mion))[0]

# marker at the end
xp = np.array(1)
yp = np.array(1)
zp = np.array(1)
line1p = axp.plot3D(xp,yp,zp,color='b',marker='o',markersize=3,lw=0)[0]

# trajectory line
line2 = axp.plot3D(xt,yt,zt,color='r',lw=1,label='ion '+"{:.2f}".format(mion2))[0]

# marker at the end
xp2 = np.array(1)
yp2 = np.array(1)
zp2 = np.array(1)
line2p = axp.plot3D(xp2,yp2,zp2,color='r',marker='o',markersize=3,lw=0)[0]

xrod = [-r0/1e-3,-r0/1e-3]
yrod = [-r0/1e-3,-r0/1e-3]
zrod = [-sa,sa]
line1rod = axp.plot3D(xrod,yrod,zrod,color='orange',lw=4,label='QMS rod')[0]
xrod1 = [-r0/1e-3,-r0/1e-3]
yrod1 = [r0/1e-3,r0/1e-3]
zrod1 = [-sa,sa]
line2rod = axp.plot3D(xrod1,yrod1,zrod1,color='orange',lw=4)[0]
xrod2 = [r0/1e-3,r0/1e-3]
yrod2 = [r0/1e-3,r0/1e-3]
zrod2 = [-sa,sa]
line3rod = axp.plot3D(xrod2,yrod2,zrod2,color='orange',lw=4)[0]
xrod3 = [r0/1e-3,r0/1e-3]
yrod3 = [-r0/1e-3,-r0/1e-3]
zrod3 = [-sa,sa]
line4rod = axp.plot3D(xrod3,yrod3,zrod3,color='orange',lw=4)[0]

# define field arrows in the plot in mm
# quiver can only generate arrows with different length by
# using individual axp.quiver commands
xgrid, ygrid, zgrid = np.meshgrid(np.arange(-sa,sa+2*sa/(gxy),2*sa/(gxy)),
                                  np.arange(-sa,sa+2*sa/(gxy),2*sa/(gxy)),
                                  np.arange(-sa*hvsw,hvsw*sa+2*hvsw*sa/gz,2*hvsw*sa/gz))

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
            uefeld[i,j,m], vefeld[i,j,m], wefeld[i,j,m] = Efield(phi0,[xb,yb,zb])
            lengthe[i,j,m] = np.sqrt(Efield(phi0,[xb,yb,zb])[0]**2+Efield(phi0,[xb,yb,zb])[1]**2+Efield(phi0,[xb,yb,zb])[2]**2)
            lengthe[i,j,m] = 0.0001
            axp.quiver(xgrid[i,j,m], ygrid[i,j,m], zgrid[i,j,m], 
                       uefeld[i,j,m], vefeld[i,j,m], wefeld[i,j,m], 
                       pivot= 'tail', color="g", length=lengthe[i,j,m]*scalee, alpha = 0.5)
            
# plotting the arrows at 0,0,0 twice to generate legend            
axp.quiver(xgrid[i,j,m], ygrid[0,0,0], zgrid[0,0,0], 
                       uefeld[0,0,0], vefeld[0,0,0], wefeld[0,0,0], 
                       pivot = 'tail', length=lengthe[0,0,0]*scalee,color="g",label="E field", alpha = 0.5)
axp.legend(fontsize=6)

# info box
infobox = ''
if model == "QMS":
    infobox += 'E0_ions: ' + "{:.2f}".format(E0) + ' (eV)\n' 
    infobox += 'm_ions: ' + "{:.2f}".format(mion) + ' (amu)\n' 
    infobox += 'U: ' + "{:.2f}".format(U) + ' (V)\n' 
    infobox += 'V: ' + "{:.2f}".format(V) + ' (V)\n'     
infobox += 'f: ' + "{:.0e}".format(f) + ' (1/s)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
axp.text2D(0,1,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=axp.transAxes)

tscan = []
vx = []
vy = []
vx2 = []
vy2 = []
ax2 = fig.add_subplot(1,2,2)
linevx = ax2.plot(tscan,vx,color='r',lw=0.5,label='vx '+"{:.2f}".format(mion))[0]
linevy = ax2.plot(tscan,vy,color='g',lw=0.5,label='vy '+"{:.2f}".format(mion))[0]
linevx2 = ax2.plot(tscan,vx2,color='orange',lw=0.5,label='vx '+"{:.2f}".format(mion2))[0]
linevy2 = ax2.plot(tscan,vy2,color='teal',lw=0.5,label='vy '+"{:.2f}".format(mion2))[0]
ax2.set(xlabel="t (microseconds)",ylabel="v_x, v_y (km/s)")
ax2.set(xlim=(0,tend/1e-6),ylim=(-vmax/1e3,vmax/1e3))
ax2.legend(fontsize=6)

fig.tight_layout(pad=3)   

# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
    global t,v,x,a,v2,x2
    global xt,yt,zt,xt2,yt2,zt2
    global vx,vy,tscan,vx2,vy2
   
    print('t '+"{:.2e}".format(t)+' s'+' of '+"{:.2e}".format(tend)+' s')    
    
    # Derivative of x and v at time t
    def deriv1(x,t):
        r = np.array([x[0],x[1],x[2]])
        vd = np.array([x[3],x[4],x[5]])
        phit = phi(U,V,t)
        accel = el/mp*Efield(phit,r)
        return np.array([vd[0],vd[1],vd[2],accel[0],accel[1],accel[2]])         
    
    ## Adding 
    time = 0
    while time<dtreplot:
        time += dt
        phit = phi(U,V,t)
        
              
        # Euler Step Particle 1
        #t += dt
        #v += el/mp * Efield(phit,x) * dt
        #x += dt * v 
        
        #rk4 step    
        xrk = [x[0],x[1],x[2],v[0],v[1],v[2]]
        half_tau = 0.5*dt
        F1 = deriv1(xrk,t)
        t_half = t + half_tau
        xtemp = xrk + half_tau*F1
        F2 = deriv1(xtemp,t_half)
        xtemp = xrk + half_tau*F2
        F3 = deriv1(xtemp,t_half)
        xtemp = xrk + half_tau*F3
        t_full = t + dt
        xtemp = xrk + dt*F3
        F4 = deriv1(xtemp,t_full)
        xrkout = xrk+dt/6*(F1+F4+2*(F2+F3))        
        x = [xrkout[0],xrkout[1],xrkout[2]]
        v = [xrkout[3],xrkout[4],xrkout[5]]
        t += dt
        
        # Euler Step Particle 2
        v2 += el/mp2 * Efield(phit,x2) * dt
        x2 += dt * v2      
                
    
    # Update the new positions in vector
    xt = np.append(xt,x[0])
    yt = np.append(yt,x[1])
    zt = np.append(zt,x[2])
    
    xt2 = np.append(xt2,x2[0])
    yt2 = np.append(yt2,x2[1])
    zt2 = np.append(zt2,x2[2])
    
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
    
    # real points in m, xgrid points in mm
    line2.set_xdata(xt2/1e-3)
    line2.set_ydata(yt2/1e-3) 
    line2.set_3d_properties(zt2/1e-3) 
    xp2 = np.array([x2[0]])
    yp2 = np.array([x2[1]])
    zp2 = np.array([x2[2]])
    line2p.set_xdata(xp2/1e-3)
    line2p.set_ydata(yp2/1e-3) 
    line2p.set_3d_properties(zp2/1e-3) 
    
    
    # rotate the view
    axp.view_init(elev,45+t/trotation*360)

    vx = np.append(vx,v[0]/1e3)
    vy = np.append(vy,v[1]/1e3)
    tscan = np.append(tscan,t/1e-6)
    
    linevx.set_xdata(tscan)
    linevx.set_ydata(vx) 
    linevy.set_xdata(tscan)
    linevy.set_ydata(vy) 
    
    vx2 = np.append(vx2,v2[0]/1e3)
    vy2 = np.append(vy2,v2[1]/1e3)    
    linevx2.set_xdata(tscan)
    linevx2.set_ydata(vx2) 
    linevy2.set_xdata(tscan)
    linevy2.set_ydata(vy2) 
    
    # Update the title of the plot
    fig.suptitle('time: ' + "{:.0f}".format(t/1e-6)+' microseconds')
  
         
 
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=180)

