#
# Geo- vs. heliozentrisches Weltbild
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

model = "Kopernikus"


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

# --------------------------------------------------------------------------
# Physical parameter Simulation
# --------------------------------------------------------------------------
if model == "Kopernikus":
  m = 1    
  animationname = 'kopernikus.mp4'


# --------------------------------------------------------------------------
# Setup Plot
# --------------------------------------------------------------------------
fig, axp = plt.subplots(1,2,figsize=(10,4.5))

axp[0].set(xlabel="x (AU)",ylabel="y (AU)")
axp[1].set(xlabel="x (AU)",ylabel="y (AU)")
scale = 4
axp[0].set(xlim=(-scale,scale),ylim=(-scale,scale))
axp[1].set(xlim=(-scale,scale),ylim=(-scale,scale))

# heliozentric
xs_h = []
ys_h = []

xe_h = []
ye_h = []

xm_h = []
ym_h = []


# trajectory lines 
lines_h = axp[0].plot(xs_h,ys_h,color='m',lw=1)[0]
linee_h = axp[0].plot(xe_h,ye_h,color='b',lw=1)[0]
linem_h = axp[0].plot(xe_h,ye_h,color='orange',lw=1)[0]

# dots
linesp_h = axp[0].plot(xs_h,ys_h,color='m',lw=0,marker='o',markersize=10,label='sun')[0]
lineep_h = axp[0].plot(xe_h,ye_h,color='b',lw=0,marker='o',markersize=6,label='earth')[0]
linemp_h = axp[0].plot(xm_h,ym_h,color='orange',lw=0,marker='o',markersize=4,label='mars')[0]

# connection line
linek = axp[0].plot(xe_h,ye_h,color='r',linestyle='dashed',lw=1)[0]

axp[0].legend(fontsize=6,loc=1)

# heliozentric
xs_g = []
ys_g = []

xe_g = []
ye_g = []

xm_g = []
ym_g = []



# trajectory lines 
lines_g = axp[1].plot(xs_h,ys_h,color='m',lw=1)[0]
linee_g = axp[1].plot(xe_h,ye_h,color='b',lw=1)[0]
linem_g = axp[1].plot(xe_h,ye_h,color='orange',lw=1)[0]

# dots
linesp_g = axp[1].plot(xs_h,ys_h,color='m',lw=0,marker='o',markersize=10,label='sun')[0]
lineep_g = axp[1].plot(xe_h,ye_h,color='b',lw=0,marker='o',markersize=6,label='earth')[0]
linemp_g = axp[1].plot(xm_h,ym_h,color='orange',lw=0,marker='o',markersize=4,label='mars')[0]

# connection line
lineg = axp[1].plot(xe_h,ye_h,color='r',linestyle='dashed',lw=1)[0]

axp[1].legend(fontsize=6,loc=1)


# info box
infobox = ''
infobox += 'heliocentric world' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
axp[0].text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=axp[0].transAxes)


# info box
infobox = ''
infobox += 'geocentric world' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
axp[1].text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=axp[1].transAxes)


   
   
# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
    global t,v,x,a
    global xt,yt,zt
   
    print("{:.3f}".format(t),' of ',"{:.3f}".format(steps*dt))
    t += dt
    
    # ------------------------------------------------------
    # heliozentrisches Weltbild
    # ------------------------------------------------------
    ex = re * np.cos(omegae*t)    
    ey = re * np.sin(omegae*t)    
   
    mx = rm * np.cos(omegam*t)    
    my = rm * np.sin(omegam*t)    
 
    sx = 0
    sy = 0
    
    linekx = [ex,mx,sx]
    lineky = [ey,my,sy]
    
    # Update the new positions in vector
    xe_h.append(ex)
    ye_h.append(ey)
    
    xm_h.append(mx)
    ym_h.append(my)
    
    xs_h.append(sx)
    ys_h.append(sy)
    
    # real points in m, xgrid points in mm
    linee_h.set_xdata(xe_h)
    linee_h.set_ydata(ye_h) 
    lineep_h.set_xdata(ex)
    lineep_h.set_ydata(ey) 
   
    linem_h.set_xdata(xm_h)
    linem_h.set_ydata(ym_h) 
    linemp_h.set_xdata(mx)
    linemp_h.set_ydata(my) 
    
    lines_h.set_xdata(xs_h)
    lines_h.set_ydata(ys_h) 
    linesp_h.set_xdata(sx)
    linesp_h.set_ydata(sy)   
   
    linek.set_xdata(linekx)
    linek.set_ydata(lineky) 
 
    
    # ------------------------------------------------------
    # geozentrisches Weltbild
    # ------------------------------------------------------
    ex = re * np.cos(omegae*t)    
    ey = re * np.sin(omegae*t)    
   
    mx = rm * np.cos(omegam*t)    
    my = rm * np.sin(omegam*t)    
 
    sx = 0
    sy = 0
    
    sx = -ex
    sy = -ey
    
    mx = mx-ex
    my = my-ey
    
    ex = 0
    ey = 0
    
    linekx = [ex,mx,sx]
    lineky = [ey,my,sy]
    
    # Update the new positions in vector
    xe_g.append(ex)
    ye_g.append(ey)
    
    xm_g.append(mx)
    ym_g.append(my)
    
    xs_g.append(sx)
    ys_g.append(sy)
    
    # real points in m, xgrid points in mm
    linee_g.set_xdata(xe_g)
    linee_g.set_ydata(ye_g) 
    lineep_g.set_xdata(ex)
    lineep_g.set_ydata(ey) 
   
    linem_g.set_xdata(xm_g)
    linem_g.set_ydata(ym_g) 
    linemp_g.set_xdata(mx)
    linemp_g.set_ydata(my) 
    
    lines_g.set_xdata(xs_g)
    lines_g.set_ydata(ys_g) 
    linesp_g.set_xdata(sx)
    linesp_g.set_ydata(sy)   
   
    lineg.set_xdata(linekx)
    lineg.set_ydata(lineky) 
 
    
   
    # Update the title of the plot
    fig.suptitle('time: ' + "{:.2f}".format(t)+' (year)')
  
         
       
# ----------------------------------------
# Create the complete time sequence
# teh total time span is tend/dtplot
# returns the rate in
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)

