#
# Parabola trajectory
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as const
from scipy.special import hyp2f1

mp = 1 # proton mass

# --------------------------------------------------------------------------
# Choose Model 
# --------------------------------------------------------------------------

model = "Parabel"


# --------------------------------------------------------------------------
# Timing Simulation
# --------------------------------------------------------------------------
t = 0
dt = 0.01 # time step
tend = 5
steps = int(tend/dt)

# --------------------------------------------------------------------------
# Physical parameter Simulation
# --------------------------------------------------------------------------
if model == "Parabel":
    animationname='parabel.mp4'

vstart = 40
angle = np.pi/6
v0x = np.cos(angle)*vstart
v0y = np.sin(angle)*vstart 
x0x = 0
x0y = 50

x2x = 100
x2y = x0y+np.tan(angle)*(x2x-x0x)
tcoll = (x2x-x0x)/v0x


# Initialize Electron Position, Velocity, Acceleration
v = np.array([v0x,v0y,0])
x = np.array([0,x0y,0])
xt = np.array(0)
yt = np.array(x0y)

v1 = np.array([v0x,v0y,0])
x1 = np.array([0,x0y,0])
xt1 = np.array(0)
yt1 = np.array(x0y)

xt1f = np.array(0)
yt1f = np.array(x0y)

v2 = np.array([0,0,0])
x2 = np.array([x2x,x2y,0])
xt2 = np.array(x2x)
yt2 = np.array(x2y)



# --------------------------------------------------------------------------
# Setup Plot
# --------------------------------------------------------------------------
fig, axp = plt.subplots(1,1,figsize=(6,6))
axp.set(xlabel="x (m)",ylabel="y (m)")
axp.set(xlim=(0,120),ylim=(0,120))

# trajectory line
line1 = axp.plot(xt,yt,color='r',lw=2,label='parabola = r(g_1) + r(g_2)')[0]
line1nog = axp.plot(xt1,yt1,color='g',lw=1,label='g_1 = 0',linestyle='dashed')[0]
line1freefall = axp.plot(xt2,yt2,color='orange',lw=1,label='g_2 = -9.81 m/s^2',linestyle='dashed')[0]
line2freefall = axp.plot(xt2,yt2,color='orange',lw=1,label='falling sphere',linestyle='dashed')[0]

# marker at the end
xp = np.array(1)
yp = np.array(1)
zp = np.array(1)
line1p = axp.plot(xp,yp,color='r',marker='o',markersize=5,lw=0)[0]
line2p = axp.plot(xp,yp,color='orange',marker='o',markersize=5,lw=0)[0]


axp.legend(fontsize=10,loc=3)

# info box
infobox = ''
infobox += 'v_start: '  + "{:.0f}".format(vstart) + ' (m/s)\n' 
infobox += 'angle: '  + "{:.0f}".format(angle*180/np.pi) + ' (deg)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
axp.text(0.05,0.95,infobox, fontsize=10,bbox=props,verticalalignment='top',transform=axp.transAxes)

def acceleration():
    return np.array([0 , -9.81 , 0])
   
# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
    global t
    global v,x
    global v1,x1
    global v2,x2
    global xt,yt
    global xt1,yt1
    global xt2,yt2
   
    #print(t)        
    t += dt
    # ------------------------------------------------
    # Boris-Bunemann Pusher
    # ------------------------------------------------
        
    ## Adding 1st half 
    v = v + acceleration() * dt 
    v1 = v1  
    v2 = v2 + acceleration() * dt 
    
    # Updating position
    x = x + dt * v      
    x1 = x1 + dt * v1      
    x2 = x2 + dt * v2      
        
    # Update the new positions in vector
    if t<tcoll:
        xt = np.append(xt,x[0])
        yt = np.append(yt,x[1])
        xt1 = np.append(xt1,x1[0])
        yt1 = np.append(yt1,x1[1])
        xt2 = np.append(xt2,x2[0])
        yt2 = np.append(yt2,x2[1])
        
        line1.set_xdata(xt)
        line1.set_ydata(yt) 
        line1nog.set_xdata(xt1)
        line1nog.set_ydata(yt1) 
        line1freefall.set_xdata([x[0],x[0]])
        line1freefall.set_ydata([x1[1],x[1]]) 
    
        line2freefall.set_xdata(xt2)
        line2freefall.set_ydata(yt2) 
        
        line1p.set_xdata(x[0])
        line1p.set_ydata(x[1]) 

        line2p.set_xdata(x2[0])
        line2p.set_ydata(x2[1]) 
    
    else:    
        xt = np.append(xt,x2[0])
        yt = np.append(yt,x2[1])
        xt1 = np.append(xt1,x1[0])
        yt1 = np.append(yt1,x1[1])
        xt2 = np.append(xt2,x[0])
        yt2 = np.append(yt2,x[1])
    
        # real points in m, xgrid points in mm
        line1.set_xdata(xt)
        line1.set_ydata(yt) 
        line1nog.set_xdata(xt1)
        line1nog.set_ydata(yt1) 
        #line1freefall.set_xdata([x[0],x[0]])
        #line1freefall.set_ydata([x1[1],x[1]]) 
    
        line2freefall.set_xdata(xt2)
        line2freefall.set_ydata(yt2) 
        
        line1p.set_xdata(x2[0])
        line1p.set_ydata(x2[1]) 

        line2p.set_xdata(x[0])
        line2p.set_ydata(x[1]) 
    
    
    # Update the title of the plot
    fig.suptitle('time: ' + "{:.2f}".format(t)+' s')
  
         
       
# ----------------------------------------
# Create the complete time sequence
# teh total time span is tend/dtplot
# returns the rate in
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)

