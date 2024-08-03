#
# Animation of an etalon
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
from matplotlib.ticker import (
    AutoLocator, AutoMinorLocator)

def delta(lamb,d,n,angle):
    return 2*np.pi/lamb*n*d*np.cos(angle)

# -----------------------------------------------------------------
# Parameter
# -----------------------------------------------------------------
R = 0.7
d = 4000*1e-9 
n = 1.5
angle = np.pi/6
lambscan = np.linspace(300*1e-9,800*1e-9,10000)
nuscan = 3e8/lambscan
deltascan = delta(lambscan,d,n,angle)

def lamb2nu(x):
    return 3e8 / x

def nu2lamb(x):
    return 3e8 / (x*1e12) / 1e-9


F = 4*R/(1-R)**2
fscan = 1/(1+F*np.sin(deltascan)**2)

fig = plt.figure(figsize=(8,4))
gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1]) 
ax = fig.add_subplot(gs[0])

line1 = ax.plot(nuscan/1e12,fscan)[0]
ax.plot([nuscan[700]/1e12,nuscan[700]/1e12],[0,1],color='red',linestyle='dashed')
ax.set_ylim(0,1)
ax.set_ylabel('T')
ax.set_xlabel(r'$\nu$' + ' (THz)')


secax = ax.secondary_xaxis('top', functions=(nu2lamb,lamb2nu))
secax.xaxis.set_minor_locator(AutoMinorLocator())
secax.set_xlabel(r'$\lambda$' +  ' (nm)')

# setup info box
infobox = ''
infobox += 'd: ' + "{:.0f}".format(d/1e-9) + ' (nm)\n' 
infobox += 'n: ' + "{:.0f}".format(n) + '\n' 
infobox += 'R: ' + "{:.2f}".format(R) + ' ' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=1) 
ax.text(0.05,0.95,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax.transAxes)

ax2 = fig.add_subplot(gs[1])
ax2.set_xlim(-2,5)
ax2.set_ylim(-4,10)
ax2.axis('off')
text1 = ax2.text(0.05,0.9,'angle',transform=ax2.transAxes,backgroundcolor='white',fontsize=8)

def betaangle(alpha,n):
   return np.arcsin(np.sin(alpha)*1/n)


dfilm = 2
filmtop = ax2.plot([],[],color='red')[0]
filmbottom = ax2.plot([],[],color='red')[0]
film = ax2.fill_between([],[],[],color='lightblue')
pattern = ax2.plot([],[],color='blue')[0]
patternfinal = ax2.plot([],[],color='blue',linestyle='dashed')[0]
linetop1 = ax2.plot([],[],color='lightblue')[0]  
linetop2 = ax2.plot([],[],color='lightblue')[0]  
linetop3 = ax2.plot([],[],color='lightblue')[0]  
linetop4 = ax2.plot([],[],color='lightblue')[0]  
linebottom1 = ax2.plot([],[],color='lightblue')[0]  
linebottom2 = ax2.plot([],[],color='lightblue')[0]  
linebottom3 = ax2.plot([],[],color='lightblue')[0]  
linebottom4 = ax2.plot([],[],color='lightblue')[0]  



def createpattern(alpha):

  xfilmbottom = [-5,0,5]
  yfilmbottom = [-5*np.tan(alpha),0,5*np.tan(alpha)]
  
  xfilmtop = [-5,0,5]
  yfilmtop = [dfilm/np.cos(alpha)-5*np.tan(alpha),0+dfilm/np.cos(alpha),dfilm/np.cos(alpha)+5*np.tan(alpha)]

  def reflbottom(x0,y0):
    beta = betaangle(alpha,n)
    dist = dfilm/np.cos(beta)
    x1 = x0-dist*np.sin(alpha-beta)
    y1 = y0+dist*np.cos(alpha-beta)
    return x1,y1

  def refltop(x0,y0):
    beta = betaangle(alpha,n)
    dist = dfilm/np.cos(beta)
    x1 = x0+dist*np.sin(alpha+beta)
    y1 = y0-dist*np.cos(alpha+beta)
    return x1,y1

  xray = [0,0]
  yray = [-5,0]
  xrayfinal = []
  yrayfinal = []
    
  linetop1x = []
  linetop1y = []
  linetop2x = []
  linetop2y = []
  linetop3x = []
  linetop3y = []

  linebottom1x = []
  linebottom1y = []
  linebottom2x = []
  linebottom2y = []
  linebottom3x = []
  linebottom3y = []
  
  linebottom3x.append(0)
  linebottom3y.append(0)
  linebottom3x.append(0+20*np.sin(2*alpha))
  linebottom3y.append(0-20*np.cos(2*alpha))
    
  xn, yn = reflbottom(xray[1],yray[1])
  xray.append(xn)
  yray.append(yn)
  linetop1x.append(xn)
  linetop1y.append(yn)
  linetop1x.append(xn)
  linetop1y.append(yn+10)
  xn, yn = refltop(xray[2],yray[2])
  xray.append(xn)
  yray.append(yn)
  linebottom1x.append(xn)
  linebottom1y.append(yn)
  linebottom1x.append(xn+20*np.sin(2*alpha))
  linebottom1y.append(yn-20*np.cos(2*alpha))
  xn, yn = reflbottom(xray[3],yray[3])
  xray.append(xn)
  yray.append(yn)
  linetop2x.append(xn)
  linetop2y.append(yn)
  linetop2x.append(xn)
  linetop2y.append(yn+10)
  xn, yn = refltop(xray[4],yray[4])
  xray.append(xn)
  yray.append(yn)
  linebottom2x.append(xn)
  linebottom2y.append(yn)
  linebottom2x.append(xn+20*np.sin(2*alpha))
  linebottom2y.append(yn-20*np.cos(2*alpha))
  xn, yn = reflbottom(xray[5],yray[5])
  xray.append(xn)
  yray.append(yn)
  linetop3x.append(xn)
  linetop3y.append(yn)
  linetop3x.append(xn)
  linetop3y.append(yn+10)
  xrayfinal.append(xn)
  yrayfinal.append(yn)
  xn, yn = refltop(xn,yn)
  xrayfinal.append(xn)
  yrayfinal.append(yn)
  
  return (xfilmbottom,yfilmbottom,xfilmtop,yfilmtop,xray,yray,
          linetop1x,linetop1y,linetop2x,linetop2y,linetop3x,linetop3y,
          linebottom1x,linebottom1y,linebottom2x,linebottom2y,linebottom3x,linebottom3y,xrayfinal,yrayfinal)
  


steps = 1000

def animate(k):
    
   global film 
   print(k) 
   angle = np.pi/6+k/steps*np.pi/6 
   deltascan = delta(lambscan,d,n,angle)
   fscan = 1/(1+F*np.sin(deltascan)**2)
   
   flambda = fscan[700]
   
   line1.set_ydata(fscan)
   text1.set_text('angle: ' + "{:.2f}".format(angle/np.pi*180)+' deg')
   
   (xfilmbottom,yfilmbottom,xfilmtop,yfilmtop,xray,yray,
      linetop1x,linetop1y,linetop2x,linetop2y,linetop3x,linetop3y,
      linebottom1x,linebottom1y,linebottom2x,linebottom2y,linebottom3x,linebottom3y,
      xrayfinal,yrayfinal) = createpattern(angle)
   linetop1.set_xdata(linetop1x)
   linetop1.set_ydata(linetop1y)
   linetop1.set_alpha(flambda)
   linetop2.set_xdata(linetop2x)
   linetop2.set_ydata(linetop2y)
   linetop2.set_alpha(flambda)
   linetop3.set_xdata(linetop3x)
   linetop3.set_ydata(linetop3y)
   linetop3.set_alpha(flambda)
   
   linebottom1.set_xdata(linebottom1x)
   linebottom1.set_ydata(linebottom1y)
   linebottom1.set_alpha(1-flambda)
   linebottom2.set_xdata(linebottom2x)
   linebottom2.set_ydata(linebottom2y)
   linebottom2.set_alpha(1-flambda)
   linebottom3.set_xdata(linebottom3x)
   linebottom3.set_ydata(linebottom3y)
   linebottom3.set_alpha(1-flambda)

   filmbottom.set_xdata(xfilmbottom)
   filmbottom.set_ydata(yfilmbottom)
   filmtop.set_xdata(xfilmtop)
   filmtop.set_ydata(yfilmtop)
   pattern.set_xdata(xray)
   pattern.set_ydata(yray)
   patternfinal.set_xdata(xrayfinal)
   patternfinal.set_ydata(yrayfinal)
   
   film.remove()
   film = ax2.fill_between(xfilmbottom,yfilmtop,yfilmbottom,color='lightblue')
  
   
  
    
    
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save('airy.mp4',fps=25,dpi=300)