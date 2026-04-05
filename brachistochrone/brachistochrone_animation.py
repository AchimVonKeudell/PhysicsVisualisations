#
# Brachistochrone
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2023
#
# Code for cycloid, line and circle for two points from
# https://scipython.com/blog/the-brachistochrone-problem/
#

import numpy as np
from scipy.optimize import newton
from scipy.integrate import quad
from scipy.interpolate import interpolate
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# acceleration
g = 9.81

# starting point, end point
x1, y1 = 0, 0
x2, y2 = 1, 0.35

# y value, lower point parabola
k = 0.6

# radius rolling cylinder
radiusroll = 0.05

# output file
animationfile = 'brachistohrone_1.gif'

# number of frames generated, time step
steps = 700
dt = 2e-3


def cycloid(x2, y2, N=1000):
    # Return the path of Brachistochrone curve from (0,0) to (x2, y2).
    # The Brachistochrone curve is the path down which a bead will fall without
    # friction between two points in the least time (an arc of a cycloid).
    # It is returned as an array of N values of (x,y) between (0,0) and (x2,y2).
    # First find theta2 from (x2, y2) numerically (by Newton-Rapheson).
    def f(theta):
        return y2/x2 - (1-np.cos(theta))/(theta-np.sin(theta))
    theta2 = newton(f, np.pi/2)

    # The radius of the circle generating the cycloid.
    R = y2 / (1 - np.cos(theta2))

    theta = np.linspace(0, theta2, N)
    x = R * (theta - np.sin(theta))
    y = R * (1 - np.cos(theta))

    # The time of travel
    T = theta2 * np.sqrt(R / g)
    return x, y, T, R

def linear(x2, y2, N=10000):
    # Return the path of a straight line from (0,0) to (x2, y2)."""

    m = y2 / x2
    x = np.linspace(0, x2, N)
    y = m*x

    # The time of travel
    T = np.sqrt(2*(1+m**2)/g/m * x2)
    return x, y, T

def func(x, f, fp):
    #The integrand of the time integral to be minimized for a path f(x)."""

    return np.sqrt((1+fp(x)**2) / (2 * g * f(x)))

def circle(x2, y2, N=10000):
    # Return the path of a circular arc between (0,0) to (x2, y2).
    # The circle used is the one with a vertical tangent at (0,0).
    
    # Circle radius
    r = (x2**2 + y2**2)/2/x2

    def f(x):
        return np.sqrt(2*r*x - x**2)
    def fp(x):
        return (r-x)/f(x)

    x = np.linspace(0, x2, N)
    y = f(x)

    # Calcualte the time of travel by numerical integration.
    T = quad(func, 0, x2, args=(f, fp))[0]
    
    return x, y, T

def parabola(x2, y2, N=10000):
     
    # Mathematica for parabola formula
    # Solve[{y1 == a (x1 - x0)^2 + k, y2 == a (x2 - x0)^2 + k}, {a, x0}]
    # Simplify[%]
    
    x1 = 0
    y1 = 0
    
    a = (1/((x1 - x2)**4))*(-2* k* (x1 - x2)**2 + x2**2* y1 - 
         2*np.sqrt((x1 - x2)**4 * (k - y1)* (k - y2)) + x2**2 * y2 + x1**2*(y1 + y2) - 
         2*x1* x2* (y1 + y2)) 
    
    x0 = -((-k*(x1 - x2)**2 + x2**2*y1 + 
             np.sqrt((x1 - x2)**4*(k - y1)*(k - y2)) + x1**2* y2 - 
             x1*x2*(y1 + y2))/((x1 - x2)*(y1 - y2)))

    def f(x):
        return a*(x-x0)**2+k
    
    x = np.linspace(0, x2, N)
    y = f(x)

    return x, y

# Plot a figure comparing the four paths.
fig, ax = plt.subplots(figsize=(6,6))

xpa, ypa = parabola(x2,y2)
lineparabola = ax.plot(xpa,ypa,lw=0.5,color='green')[0]
lineparabolaroll = ax.plot(x1,y2,lw=2,color='green')[0]
lineparabolaarrow = ax.plot(x1,y2,lw=1,color='green')[0]
lineparabolarollpt = ax.plot(x1,y2,marker='o',markersize=5,color='green')[0]
lineparabolapt = ax.plot(x1,x2,marker='x',markersize=5,color='green',label='parabola')[0]
parabola = interpolate.interp1d(xpa,ypa)

xc, yc, T, R = cycloid(x2,y2)
linecycloid = ax.plot(xc,yc,lw=0.5,color='b')[0]
linecycloidroll = ax.plot(x1,y2,lw=2,color='b')[0]
linecycloidrollpt = ax.plot(x1,y2,marker='o',markersize=5,color='b')[0]
linecycloidpt = ax.plot(x1,x2,marker='x',markersize=5,color='b',label='cycloid')[0]
cycloid = interpolate.interp1d(xc,yc)

xl, yl, T = linear(x2,y2)
linelinear = ax.plot(xl,yl,lw=0.5,color='orange')[0]
linelinearroll = ax.plot(x1,y2,lw=2,color='orange')[0]
linelinearrollpt = ax.plot(x1,y2,marker='o',markersize=5,color='orange')[0]
linelinearpt = ax.plot(x1,x2,marker='x',markersize=5,label='linear',color='orange')[0]
linear = interpolate.interp1d(xl,yl)

xr, yr, T = circle(x2,y2)
linecircle = ax.plot(xr,yr,lw=0.5,color='r')[0]
linecircleroll = ax.plot(x1,y2,lw=2,color='r')[0]
linecirclerollpt = ax.plot(x1,y2,marker='o',color='r',markersize=5)[0]
linecirclept = ax.plot(x1,x2,marker='x',markersize=5,label='circle',color='r')[0]
circle = interpolate.interp1d(xr,yr)

linecycloidcirclebase = ax.plot([x1,x2],[y1,y1],lw=1,color='b')[0]
linecycloidcircle = ax.plot(xc,yc,lw=1,color='b',linestyle='dotted')[0]
linecycloidcirclept = ax.plot(xc,yc,lw=1,color='b',marker='x')[0]

ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_xlim(-0.1, 1.1)
ax.set_ylim(0.9, -0.3)
ax.set_yticks([])
ax.legend()

t = 0
microsteps = 1
deltax = 1e-5

xpap, ypap = x1, y1
vpap = 0
rollthetaparabola = np.pi
timeparabola = 0

xcp, ycp = x1, y1
vcp = 0
rollthetacycloid = np.pi
timecycloid = 0

xlp, ylp = x1, y1
vlp = 0
rollthetalinear = np.pi
timelinear = 0

xrp, yrp = x1, y1
vrp = 0
rollthetacircle = np.pi
timecircle = 0


def rollcircle(xo,yo):
    theta = np.linspace(0,2*np.pi,90)
    xroll = xo + radiusroll*np.cos(theta)
    yroll = yo + radiusroll*np.sin(theta)
    return xroll, yroll

def rollcycloid(xo,yo,R):
    theta = np.linspace(0,2*np.pi,90)
    xroll = xo + R*np.cos(theta)
    yroll = yo + R*np.sin(theta)
    return xroll, yroll

def animate(k):
    
    print(k)
    
    global xpap,ypap
    global xcp,ycp
    global xlp,ylp
    global xrp,yrp
    
    global vpap
    global vcp
    global vlp
    global vrp
    
    global rollthetaparabola
    global rollthetacircle
    global rollthetalinear
    global rollthetacycloid
    
    global timeparabola
    global timecycloid
    global timecircle
    global timelinear
    
    # parabola
    if xpap>=0 and xpap<=x2-deltax:

      xpap0 = xpap
      ypap0 = ypap 
  
      if xpap+deltax<=x2:
        dy = parabola(xpap+deltax)-parabola(xpap)
        dx = deltax
        forceangle = np.arctan(dy/dx)    
      
      apap = g*np.cos(np.pi/2-forceangle)
      vpap += apap*dt
        
      xpap += np.cos(forceangle)*vpap*dt
      ypap += np.sin(forceangle)*vpap*dt
 
        
          
      xo = xpap + np.cos(np.pi/2-forceangle)*radiusroll
      yo = ypap - np.sin(np.pi/2-forceangle)*radiusroll
      xroll, yroll = rollcircle(xo,yo)
      lineparabolaroll.set_xdata(xroll)
      lineparabolaroll.set_ydata(yroll)
  
      
      if xpap<=x2:
    
        lineparabolapt.set_xdata(xpap)
        lineparabolapt.set_ydata(ypap)
        
        rollthetaparabola += np.sqrt((xpap-xpap0)**2 + (ypap-ypap0)**2)/(2*np.pi*radiusroll) * np.pi*2
        rollxpt = xo + np.cos(rollthetaparabola)*radiusroll
        rollypt = yo + np.sin(rollthetaparabola)*radiusroll
        lineparabolarollpt.set_xdata(rollxpt)
        lineparabolarollpt.set_ydata(rollypt)
        
        timeparabola += dt*microsteps
    
    # cycloid
    if xcp>=0 and xcp<=x2-deltax:
        
      xcp0 = xcp
      ycp0 = ycp        
        
      if xcp+deltax<=x2:   
        dy = cycloid(xcp+deltax)-cycloid(xcp)
        dx = deltax
        forceangle = np.arctan(dy/dx)  
       
      acp = g*np.cos(np.pi/2-forceangle)
      vcp += acp*dt
        
      xcp += np.cos(forceangle)*vcp*dt
      ycp += np.sin(forceangle)*vcp*dt
           
      xo = xcp + np.cos(np.pi/2-forceangle)*radiusroll
      yo = ycp - np.sin(np.pi/2-forceangle)*radiusroll
      xroll, yroll = rollcircle(xo,yo)
      linecycloidroll.set_xdata(xroll)
      linecycloidroll.set_ydata(yroll)
      
      if xcp <= x2:
       
        linecycloidpt.set_xdata(xcp)
        linecycloidpt.set_ydata(ycp)
        
        rollthetacycloid += np.sqrt((xcp-xcp0)**2 + (ycp-ycp0)**2)/(2*np.pi*radiusroll) * np.pi*2
        rollxpt = xo + np.cos(rollthetacycloid)*radiusroll
        rollypt = yo + np.sin(rollthetacycloid)*radiusroll
        linecycloidrollpt.set_xdata(rollxpt)
        linecycloidrollpt.set_ydata(rollypt)
        
        shiftR = R**2-(cycloid(xcp)-R)**2
        if shiftR<0: shiftR=0
        if forceangle<=0:
          circlex, circley = rollcycloid(xcp-np.sqrt(shiftR),y1+R,R)
        else:
          circlex, circley = rollcycloid(xcp+np.sqrt(shiftR),y1+R,R)
        linecycloidcircle.set_xdata(circlex)
        linecycloidcircle.set_ydata(circley)
       
        if forceangle<=0:
          linecycloidcirclept.set_xdata([xcp,xcp-np.sqrt(shiftR)])
        else:
          linecycloidcirclept.set_xdata([xcp,xcp+np.sqrt(shiftR)])
        linecycloidcirclept.set_ydata([ycp,y1+R])
        
        timecycloid += dt*microsteps
    
    # linear
    if xlp>=0 and xlp<=x2-deltax:
     
      xlp0 = xlp
      ylp0 = ylp 
       
      if xlp+deltax<=x2:
        dy = linear(xlp+deltax)-linear(xlp)
        dx = deltax
        forceangle = np.arctan(dy/dx)    
      
      alp = g*np.cos(np.pi/2-forceangle)
      vlp += alp*dt
        
      xlp += np.cos(forceangle)*vlp*dt
      ylp += np.sin(forceangle)*vlp*dt
                  
        
      xo = xlp + np.cos(np.pi/2-forceangle)*radiusroll
      yo = ylp - np.sin(np.pi/2-forceangle)*radiusroll
      xroll, yroll = rollcircle(xo,yo)
      linelinearroll.set_xdata(xroll)
      linelinearroll.set_ydata(yroll)
  
      
      if xlp<=x2:
     
        linelinearpt.set_xdata(xlp)
        linelinearpt.set_ydata(ylp)
        
        rollthetalinear += np.sqrt((xlp-xlp0)**2 + (ylp-ylp0)**2)/(2*np.pi*radiusroll) * np.pi*2
        rollxpt = xo + np.cos(rollthetalinear)*radiusroll
        rollypt = yo + np.sin(rollthetalinear)*radiusroll
        linelinearrollpt.set_xdata(rollxpt)
        linelinearrollpt.set_ydata(rollypt)
        
        timelinear += dt*microsteps
     
   
    # circle
    if xrp>=0 and xrp<=x2-deltax:
        
      yrp0 = yrp  
      xrp0 = xrp
     
      if xrp+deltax<=x2:  
        dy = circle(xrp+deltax)-circle(xrp)
        dx = deltax
        forceangle = np.arctan(dy/dx)    
      
      arp = g*np.cos(np.pi/2-forceangle)
      vrp += arp*dt
        
      xrp += np.cos(forceangle)*vrp*dt
      yrp += np.sin(forceangle)*vrp*dt
        
          
      xo = xrp + np.cos(np.pi/2-forceangle)*radiusroll
      yo = yrp - np.sin(np.pi/2-forceangle)*radiusroll
      xroll, yroll = rollcircle(xo,yo)
      linecircleroll.set_xdata(xroll)
      linecircleroll.set_ydata(yroll)
  
      
      if xrp<=x2:
     
        linecirclept.set_xdata(xrp)
        linecirclept.set_ydata(yrp)
     
        rollthetacircle += np.sqrt((xrp-xrp0)**2 + (yrp-yrp0)**2)/(2*np.pi*radiusroll) * np.pi*2
        rollxpt = xo + np.cos(rollthetacircle)*radiusroll
        rollypt = yo + np.sin(rollthetacircle)*radiusroll
        linecirclerollpt.set_xdata(rollxpt)
        linecirclerollpt.set_ydata(rollypt)
        
        timecircle += dt*microsteps   
      
        
    fig.suptitle('t$_{cycloid}$: ' + "{:.2f}".format(timecycloid)+' s, '+
                 't$_{linear}$: ' + "{:.2f}".format(timelinear)+' s, '+
                 't$_{circle}$: ' + "{:.2f}".format(timecircle)+' s, '+
                 't$_{parabola}$: ' + "{:.2f}".format(timeparabola)+' s')
   
      

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationfile,fps=25,dpi=300)