#
# Planet orbits as a 3 body problem
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as const


# --------------------------------------------------------------------------
# Choose Model 
# --------------------------------------------------------------------------

model = "planeten3body"

# --------------------------------------------------------------------------
# Timing Simulation
# --------------------------------------------------------------------------
t = 0
dt = 1e-3 # time step
microsteps = 100
tend = 500
steps = int(tend/(microsteps*dt)) # simulation steps

Te = 1
re = 1

rm = 1.67*re
Tm = np.sqrt(Te*rm**3/re**3)

omegae = 1/Te*2*np.pi
omegam = 1/Tm*2*np.pi

G = 6.67e-11 # gravitational constant

m1 = 5.972e27 # ,ass earth
m2 = 3.98847e27 # mass sun
m3 = 1.98847e28 # mass sun

au = 1.4959e11 # astronomical unit
v1 = 2978 # m/s 
v2 = 2650 # m/s
v3 = 2650 # m/s

sproa = 3.154e7 # seconds per year

def accel(x1,y1,x2,y2,x3,y3):
    # planet 1
    # force to body 2
    r2 = ((x1-x2)*au)**2+((y1-y2)*au)**2
    r = np.sqrt(r2)
    acc = -G*m2/r2    
    ax1 = (x1-x2)*au/r*acc
    ay1 = (y1-y2)*au/r*acc
    # force to body 3
    r2 = ((x1-x3)*au)**2+((y1-y3)*au)**2
    r = np.sqrt(r2)
    acc = -G*m3/r2    
    ax1 += (x1-x3)*au/r*acc
    ay1 += (y1-y3)*au/r*acc

    # planet 2
    # force to body 1
    r2 = ((x2-x1)*au)**2+((y2-y1)*au)**2
    r = np.sqrt(r2)
    acc = -G*m1/r2    
    ax2 = (x2-x1)*au/r*acc
    ay2 = (y2-y1)*au/r*acc
    # force to body 3
    r2 = ((x2-x3)*au)**2+((y2-y3)*au)**2
    r = np.sqrt(r2)
    acc = -G*m3/r2    
    ax2 += (x2-x3)*au/r*acc
    ay2 += (y2-y3)*au/r*acc

    # planet 3
    # force to body 1
    r2 = ((x3-x1)*au)**2+((y3-y1)*au)**2
    r = np.sqrt(r2)
    acc = -G*m1/r2    
    ax3 = (x3-x1)*au/r*acc
    ay3 = (y3-y1)*au/r*acc
    # force to body 2
    r2 = ((x3-x2)*au)**2+((y3-y2)*au)**2
    r = np.sqrt(r2)
    acc = -G*m2/r2    
    ax3 += (x3-x2)*au/r*acc
    ay3 += (y3-y2)*au/r*acc
    
    return ax1, ay1, ax2, ay2, ax3, ay3
    
# --------------------------------------------------------------------------
# Physical parameter Simulation
# --------------------------------------------------------------------------
if model == "planeten3body":
  m = 1    
  animationname = 'planets3body.mp4'


# --------------------------------------------------------------------------
# Setup Plot
# --------------------------------------------------------------------------
fig, axp = plt.subplots(1,1,figsize=(4.5,4.5))

axp.set(xlabel="x (AU)",ylabel="y (AU)")
scale = 10
axp.set(xlim=(-scale,scale),ylim=(-scale,scale))

# trajectories
x1_h = []
y1_h = []

x2_h = []
y2_h = []

x3_h = []
y3_h = []


# trajectory lines 
line1_h = axp.plot(x1_h,y1_h,color='m',lw=0.5,alpha=0.7)[0]
line2_h = axp.plot(x2_h,y2_h,color='b',lw=0.5,alpha=0.7)[0]
line3_h = axp.plot(x3_h,y3_h,color='orange',lw=0.5,alpha=0.7)[0]

# dots
axp.plot(0,0,color='r',lw=0,marker='+',markersize=6,label='center of mass')
line1p_h = axp.plot(x1_h,y1_h,color='m',lw=0,marker='o',markersize=m1/3e27,label='planet 1')[0]
line2p_h = axp.plot(x2_h,y2_h,color='b',lw=0,marker='o',markersize=m2/3e27,label='planet 2')[0]
line3p_h = axp.plot(x3_h,y3_h,color='orange',lw=0,marker='o',markersize=m3/3e27,label='planet 3')[0]

quiver1 = axp.quiver(0,0, 0, 1, pivot = 'tail', color="b", alpha=0.8)
quivera1 = axp.quiver(0,0, 0, 1, pivot = 'tail', color="r", alpha=0.8)

quiver2 = axp.quiver(0,0, 0, 1, pivot = 'tail', color="b", alpha=0.8)
quivera2 = axp.quiver(0,0, 0, 1, pivot = 'tail', color="r", alpha=0.8)

quiver3 = axp.quiver(0,0, 0, 1, pivot = 'tail', color="b", alpha=0.8)
quivera3 = axp.quiver(0,0, 0, 1, pivot = 'tail', color="r", alpha=0.8)

axp.legend(fontsize=8,loc=1)


# info box
#infobox = ''
#infobox += 'heliocentric world' 
#props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
#axp.text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=axp.transAxes)

# initial conditions
# position in au units
# velocity in real units
# time in years
x1 = 0.2
y1 = 1 

vy1 = v1*0.2 
vx1 = 0

x2 = 1.4
y2 = 0 

vy2 = v1*0.6 
vx2 = 0

x3 = -1.2
y3 = 0 

vy3 = -v1*0.8 
vx3 = 0

vsx = (vx1*m1+vx2*m2+vx3*m3)/(m1+m2+m3)
vsy = (vy1*m1+vy2*m2+vy3*m3)/(m1+m2+m3)

sx = (x1*m1+x2*m2+x3*m3)/(m1+m2+m3)
sy = (y1*m1+y2*m2+y3*m3)/(m1+m2+m3)

x1 = x1 - sx
x2 = x2 - sx
x3 = x3 - sx
y1 = y1 - sy
y2 = y2 - sy
y3 = y3 - sy
vx1 = vx1 - vsx
vx2 = vx2 - vsx
vx3 = vx3 - vsx
vy1 = vy1 - vsy
vy2 = vy2 - vsy
vy3 = vy3 - vsy
   
# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
    
    global t,v,x,a
    global xt,yt,zt
    global x1,y1,vx1,vy1
    global x2,y2,vx2,vy2
    global x3,y3,vx3,vy3
    global quiver1, quivera1
    global quiver2, quivera2
    global quiver3, quivera3
    
    print("{:.3f}".format(t),' of ',"{:.3f}".format(steps*dt*microsteps))
    t += dt*microsteps
    
    
    # ------------------------------------------------------
    # planeten
    # ------------------------------------------------------
    for i in range(microsteps):
      a1x, a1y, a2x, a2y, a3x, a3y = accel(x1,y1,x2,y2,x3,y3)
    
      vx1 = vx1 + a1x*dt*sproa
      vy1 = vy1 + a1y*dt*sproa
    
      x1 = x1 + vx1/au*dt*sproa    
      y1 = y1 + vy1/au*dt*sproa
   
      vx2 = vx2 + a2x*dt*sproa
      vy2 = vy2 + a2y*dt*sproa
    
      x2 = x2 + vx2/au*dt*sproa    
      y2 = y2 + vy2/au*dt*sproa
   
      vx3 = vx3 + a3x*dt*sproa
      vy3 = vy3 + a3y*dt*sproa
    
      x3 = x3 + vx3/au*dt*sproa    
      y3 = y3 + vy3/au*dt*sproa
   
  
    # Update the new positions in vector
    x1_h.append(x1)
    y1_h.append(y1)
    
    #print(a1x,a1y)
    
    line1_h.set_xdata(x1_h)
    line1_h.set_ydata(y1_h) 
    line1p_h.set_xdata(x1)
    line1p_h.set_ydata(y1) 
   
    x2_h.append(x2)
    y2_h.append(y2)
    
    line2_h.set_xdata(x2_h)
    line2_h.set_ydata(y2_h) 
    line2p_h.set_xdata(x2)
    line2p_h.set_ydata(y2) 
   
    x3_h.append(x3)
    y3_h.append(y3)
    
    line3_h.set_xdata(x3_h)
    line3_h.set_ydata(y3_h) 
    line3p_h.set_xdata(x3)
    line3p_h.set_ydata(y3) 
   
    

    # ------------------------------------------------------
    # arrows velocity, acceleration
    # ------------------------------------------------------    
    
    a1x, a1y, a2x, a2y, a3x, a3y = accel(x1,y1,x2,y2,x3,y3)
    quiver1.remove()
    quiver1 = axp.quiver(x1,y1, vx1, vy1, scale = 1e5, pivot = 'tail', color="b", alpha=0.8, label="velocity")
    quivera1.remove()
    quivera1 = axp.quiver(x1,y1, a1x, a1y, scale = 1e-2, pivot = 'tail', color="r", alpha=0.8, label="acceleration")
    quiver2.remove()
    quiver2 = axp.quiver(x2,y2, vx2, vy2, scale = 1e5, pivot = 'tail', color="b", alpha=0.8, label="velocity")
    quivera2.remove()
    quivera2 = axp.quiver(x2,y2, a2x, a2y, scale = 1e-2, pivot = 'tail', color="r", alpha=0.8, label="acceleration")
    quiver3.remove()
    quiver3 = axp.quiver(x3,y3, vx3, vy3, scale = 1e5, pivot = 'tail', color="b", alpha=0.8, label="velocity")
    quivera3.remove()
    quivera3 = axp.quiver(x3,y3, a3x, a3y, scale = 1e-2, pivot = 'tail', color="r", alpha=0.8, label="acceleration")

        
   
    # Update the title of the plot
    fig.suptitle('time: ' + "{:.2f}".format(t)+' (year)')
  
         
       
# ----------------------------------------
# Create the complete time sequence
# teh total time span is tend/dtplot
# returns the rate in
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)

