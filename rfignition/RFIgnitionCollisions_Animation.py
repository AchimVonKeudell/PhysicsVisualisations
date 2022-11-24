#
# Ignition in a RF discharge
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

el = 1.6e-19
me = 9e-31

dt = 1e-9
omega = 13.56*1e6*2*np.pi
steps = 500
E0 = 100
nu = 1e8

# Initialize Velocities
vx = 0
vy = 0
vz = 0
x = 0
y = 0
ax = 0
ay = 0

vx1 = 0
vy1 = 0
vz1 = 0
x1 = 0
y1 = 0
ax1 = 0
ay1 = 0

t = 0
ener = [0]
ener1 = [0]
tscan =[0]

fig, axp = plt.subplots(1,1,figsize=(8,4.5))
line1 = axp.plot(tscan,ener,color='r',lw=1, label ='no collisions')[0]
line2 = axp.plot(tscan,ener1,color='b',lw=1, label ='with collisions')[0]
axp.set(xlabel="time (microseconds)",ylabel="energy (eV)")
axp.set(xlim=(0,steps*10*dt/1e-6),ylim=(0,50))
axp.legend(fontsize=8,labelspacing=0.1,facecolor='lightblue')
   
def animate(k):
    global t, vx,vy,vz,x,y,ax,ay,vx1,vy1,vz1,ax1,ay1,x1,y1
    for i in range(10):
        t += dt
        tscan.append(t/1e-6)
        Ex = np.sin(omega*t)*E0
        Ey = 0

        # ------------------------------------------------
        # First particle no collisions
        # ------------------------------------------------
        # 1st Half Kick
        vx += ax*0.5*dt
        vy += ay*0.5*dt
 
        # Move Partciles 
        x += vx*dt
        y += vy*dt

        # Acceleration  
        ax = Ex*el/me
        ay = Ey*el/me
   
        # 2nd Half Kick
        vx += ax*0.5*dt;
        vy += ay*0.5*dt;
    
        ener.append(0.5*me/el*(vx**2+vy**2+vz**2))
    
 
        # ----------------------------------------------
        # Second particle collisions
        # ----------------------------------------------

        # Half Kick
        vx1 += ax1*0.5*dt
        vy1 += ay1*0.5*dt
 
        # Move Partciles 
        x1 += vx1*dt
        y1 += vy1*dt

        # Acceleration  
        ax1 = Ex*el/me
        ay1 = Ey*el/me

        # ------------------------------------------------------
        # Collision calculation
        # from eduPIC code example by Z. Donko et al.
        # ------------------------------------------------------
        if np.random.rand() < 1-np.exp(-nu*dt):
            v = np.sqrt(vx1**2+vy1**2+vz1**2)
            # calulate Euler angles
            if vx1 == 0:
                theta = 0.5*np.pi
            else:
                theta = np.arctan(np.sqrt(vy1**2+vz1**2)/vx1)
            if vy1 == 0:
                if vz1 > 0:
                    phi = 0.5*np.pi
                else:
                    phi = -0.5*np.pi
            else:        
                phi = np.arctan(vz1/vy1)

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
            vx1 = v * (ct * cc - st * sc * ce);
            vy1 = v * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
            vz1 = v * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    
        # 2nd half kick
        vx1 += ax1*0.5*dt;
        vy1 += ay1*0.5*dt;
    
        ener1.append(0.5*me/el*(vx1**2+vy1**2+vz1**2))

        line1.set_xdata(tscan)
        line1.set_ydata(ener) 
        line2.set_xdata(tscan)
        line2.set_ydata(ener1) 
         
       
# ----------------------------------------
# Create the complete time sequence
# teh total time span is tend/dtplot
# returns the rate in
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save('rfignition.gif',fps=25,dpi=180)

