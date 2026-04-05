#
# Ion-induced surface roughening
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

from numba import jit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.ndimage.filters import gaussian_filter

model = 'normalincidence'
model = 'glancingincidence'

ionrange = 1
lam = 2

steps = 1000
ministeps = 1e7
blurradius = 5
NGrid = 600
fringe  = 2

if model == 'normalincidence':
   animationfile = 'ionroughnormal.mp4'
   angle = 0
   ionrangey = ionrange*np.sin(angle)
   ionrangez = ionrange*np.cos(angle)

if model == 'glancingincidence':
   animationfile = 'ionroughglancing.mp4'
   angle = 45*np.pi/180
   ionrangey = ionrange*np.sin(angle)
   ionrangez = ionrange*np.cos(angle)

surface = np.zeros([NGrid,NGrid])
vmaximum = np.sqrt(steps*ministeps/NGrid**2)

x = np.arange(0, NGrid-2*fringe, 1)
y = np.arange(0, NGrid-2*fringe, 1)
x, y = np.meshgrid(x, y)

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, surface[fringe:NGrid-fringe,fringe:NGrid-fringe], vmin= -vmaximum, vmax = vmaximum,cmap='viridis',edgecolor='None')

#cbar = fig.colorbar(surface[fringe:NGrid-fringe,fringe:NGrid-fringe],ax=ax,label='height')
#cbar.set_ticks([])

# info box
infobox = ''
infobox += 'ion range: ' + "{:.0f}".format(ionrange) + ' (a.u.)\n' 
infobox += 'ion angle: ' + "{:.0f}".format(angle/np.pi*180) + ' (deg)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 


@jit(nopython = True)
def evolve(surface):
  
  d = []
  d.append([1,0])
  d.append([-1,0])
  d.append([0,1])
  d.append([0,-1])
  d.append([2,0])
  d.append([-2,0])
  d.append([0,2])
  d.append([0,-2])
  d.append([1,1])
  d.append([-1,-1])
  d.append([1,-1])
  d.append([1,1])
  d.append([-1,2])
  d.append([1,2])
  d.append([-1,-2])
  d.append([1,-2])
  d.append([-2,1])
  d.append([2,1])
  d.append([-2,-1])
  d.append([1,-1])
  
  NBneighbors = len(d)
   
  for k in range(int(ministeps)):
  
    xi = int(fringe+np.random.rand()*(NGrid-2*fringe))
    yi = int(fringe+np.random.rand()*(NGrid-2*fringe))
    
    x1 = xi
    y1 = yi
    z1 = surface[xi,yi]
    
    for m in range(NBneighbors):
        
      x2 = xi+d[m][0]
      y2 = yi+d[m][1]
      z2 = surface[x2,y2]
      
      de = np.sqrt((z2-z1-ionrangez)**2+(x2-x1)**2+(y2-y1+ionrangey)**2)
      dz = np.exp(-de/lam)
      surface[xi,yi] -= dz 
  
  h = np.mean(surface[fringe:NGrid-fringe,fringe:NGrid-fringe])
  surface -= h  

  
def animate(f):
  global surface
 
  print(f)
  evolve(surface)
  
  ax.clear()
  ax.set(zlim=(-vmaximum,vmaximum))
  ax.set_axis_off()
  
  blurred = gaussian_filter(surface[fringe:NGrid-fringe,fringe:NGrid-fringe], sigma=blurradius)
  #ax.quiver(NGrid/2, NGrid/2, 2*vmaximum, 0, ionrangey*NGrid/(2*vmaximum), -ionrangez, length=NGrid/10, color='red', normalize=True)
  ax.plot_surface(x, y, blurred, cmap='seismic',alpha=0.5)
  ax.view_init(elev=90*np.sin(np.pi*f/steps), azim=45+f)
  
  
  ax.text2D(0,1,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)
  fig.suptitle('Ions: ' + "{:.0f}".format(ministeps*(f+1)/1e4)+' x 10$^4$')


anim = animation.FuncAnimation(fig,animate,interval=0,frames=steps)
anim.save(animationfile,fps=25,dpi=300) 

