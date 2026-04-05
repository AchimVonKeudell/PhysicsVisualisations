#
# Planet orbits
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

model = "planeten"


# --------------------------------------------------------------------------
# Timing Simulation
# --------------------------------------------------------------------------
t = 0
dt = 1e-2 # time step
steps = 2000 # simulation steps

Te = 1
re = 1

rm = 1.67*re
Tm = np.sqrt(Te*rm**3/re**3)

omegae = 1/Te*2*np.pi
omegam = 1/Tm*2*np.pi

G = 6.67e-11 # gravitational constant
ms = 1.98847e30 # mass sun
me = 5.972e24 # ,ass earth
au = 1.4959e11 # astronomical unit
ve = 29780 # m/s 
vm = 26500 # m/s
sproa = 3.154e7 # seconds per year

def accel(x,y):
    r2 = (x*au)**2+(y*au)**2
    r = np.sqrt(r2)
    acc = -G*ms/r2    
    ax = x*au/r*acc
    ay = y*au/r*acc
    return ax, ay
    
# --------------------------------------------------------------------------
# Physical parameter Simulation
# --------------------------------------------------------------------------
if model == "planeten":
  m = 1    
  animationname = 'planets.mp4'


# --------------------------------------------------------------------------
# Setup Plot
# --------------------------------------------------------------------------
fig, axp = plt.subplots(1,1,figsize=(4.5,4.5))

axp.set(xlabel="x (AU)",ylabel="y (AU)")
scale = 6
axp.set(xlim=(-scale,scale),ylim=(-scale,scale))

# heliozentric
xs_h = []
ys_h = []

xe_h = []
ye_h = []

xm_h = []
ym_h = []


# trajectory lines 
lines_h = axp.plot(xs_h,ys_h,color='m',lw=1)[0]
linee_h = axp.plot(xe_h,ye_h,color='b',lw=1)[0]
linem_h = axp.plot(xe_h,ye_h,color='orange',lw=1)[0]

# dots
linesp_h = axp.plot(xs_h,ys_h,color='m',lw=0,marker='o',markersize=10,label='sun')[0]
lineep_h = axp.plot(xe_h,ye_h,color='b',lw=0,marker='o',markersize=6,label='planet 1')[0]
linemp_h = axp.plot(xm_h,ym_h,color='orange',lw=0,marker='o',markersize=4,label='planet 2')[0]

quivere = axp.quiver(0,0, 0, 1, pivot = 'tail', color="b", alpha=0.8, label="velocity")
quiverae = axp.quiver(0,0, 0, 1, pivot = 'tail', color="r", alpha=0.8, label="acceleration")

quiverm = axp.quiver(0,0, 0, 1, pivot = 'tail', color="b", alpha=0.8)
quiveram = axp.quiver(0,0, 0, 1, pivot = 'tail', color="r", alpha=0.8)


axp.legend(fontsize=6,loc=1)


# info box
#infobox = ''
#infobox += 'heliocentric world' 
#props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
#axp.text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=axp.transAxes)

# initial conditions
# position in au units
# velocity in real units
# time in years
ex = 1
ey = 0 

vey = ve*1.3 
vex = 0

mx = -1.67
my = 0 

vmy = -vm*1.05 
vmx = 0
   
# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
    
    global t,v,x,a
    global xt,yt,zt
    global ex,ey,vex,vey
    global quivere, quiverae
    global mx,my,vmx,vmy
    global quiverm, quiveram
    
    print("{:.3f}".format(t),' of ',"{:.3f}".format(steps*dt))
    t += dt
    
    # ------------------------------------------------------
    # earth
    # ------------------------------------------------------
    aex, aey = accel(ex,ey)
    vex = vex + aex*dt*sproa
    vey = vey + aey*dt*sproa
    
    ex = ex + vex/au*dt*sproa    
    ey = ey + vey/au*dt*sproa
   
    # Update the new positions in vector
    xe_h.append(ex)
    ye_h.append(ey)
    
    linee_h.set_xdata(xe_h)
    linee_h.set_ydata(ye_h) 
    lineep_h.set_xdata(ex)
    lineep_h.set_ydata(ey) 
   
    # ------------------------------------------------------
    # mars
    # ------------------------------------------------------
    amx, amy = accel(mx,my)
    vmx = vmx + amx*dt*sproa
    vmy = vmy + amy*dt*sproa
    
    mx = mx + vmx/au*dt*sproa    
    my = my + vmy/au*dt*sproa
   
    # Update the new positions in vector
    xm_h.append(mx)
    ym_h.append(my)
    
    linem_h.set_xdata(xm_h)
    linem_h.set_ydata(ym_h) 
    linemp_h.set_xdata(mx)
    linemp_h.set_ydata(my) 
   
    # ------------------------------------------------------
    # sun
    # ------------------------------------------------------
    sx = 0
    sy = 0
          
    xs_h.append(sx)
    ys_h.append(sy)
    
    lines_h.set_xdata(xs_h)
    lines_h.set_ydata(ys_h) 
    linesp_h.set_xdata(sx)
    linesp_h.set_ydata(sy)   
    

    # ------------------------------------------------------
    # arrows velocity, acceleration
    # ------------------------------------------------------    
    quivere.remove()
    quivere = axp.quiver(ex,ey, vex, vey, scale = 1e5, pivot = 'tail', color="b", alpha=0.8, label="velocity")

    aex, aey = accel(ex,ey)   
    quiverae.remove()
    quiverae = axp.quiver(ex,ey, aex, aey, scale = 1e-2, pivot = 'tail', color="r", alpha=0.8, label="acceleration")

    quiverm.remove()
    quiverm = axp.quiver(mx,my, vmx, vmy, scale = 1e5, pivot = 'tail', color="b", alpha=0.8, label="velocity")

    amx, amy = accel(mx,my)   
    quiveram.remove()
    quiveram = axp.quiver(mx,my, amx, amy, scale = 1e-2, pivot = 'tail', color="r", alpha=0.8, label="acceleration")
        
   
    # Update the title of the plot
    fig.suptitle('time: ' + "{:.2f}".format(t)+' (year)')
  
         
       
# ----------------------------------------
# Create the complete time sequence
# teh total time span is tend/dtplot
# returns the rate in
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)

