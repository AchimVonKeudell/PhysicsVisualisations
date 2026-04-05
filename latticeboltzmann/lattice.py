#
# Flow phenomena
# Simulation via Lattice Boltzmann Method
# Examples for physics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
# Code adapted from:
# Create Your Own Lattice Boltzmann Simulation (With Python)
# Philip Mocz (2020) Princeton Univeristy, @PMocz
#
#

#from numba import njit
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import math

#model = 'narrow'
#model = 'cylinder'
#model = 'square'
model = 'tube'
model = 'ellipse'

# Simulation parameters
Nx = 600    # resolution x-dir
Ny = 200    # resolution y-dir
rho0 = 100    # average density
obstaclex = 100
obstacley = 100

if model == 'cylinder':
  animationfile = 'latticecylinder_long_2.mp4'  
  obstacleradius = 15
  tau = 0.6    # collision timescale
  show = 'vorticity'
  cmin = -0.015
  cmax = 0.015
  #show = 'density'
  plotstream = False
  plotarrow = True
  Nt = 9000   # number of timesteps
  circlex = [100,160]
  circley = [140,60]

if model == 'ellipse':
  animationfile = 'latticeellipse_8.mp4'  
  obstacleradiusa = 70
  obstacleradiusb = 10
  tau = 0.6    # collision timescale
  show = 'vorticity'
  cmin = -0.015
  cmax = 0.015
  #show = 'density'
  show = 'velocity'
  cmin = 0
  cmax = 30
  stallangle = 15
  plotstream = False
  plotarrow = True
  Nt = 9000   # number of timesteps
  circlex = [150]
  circley = [int(Ny/2)]

if model == 'square':
  animationfile = 'latticesquare.mp4'  
  obstacleheight = 30
  obstaclewidth = 5
  tau = 0.6    # collision timescale
  show = 'vorticity'
  cmin = -0.015
  cmax = 0.015
  #show = 'density'
  #cmin = -0.01
  #cmax = 0.01  
  plotstream = False
  plotarrow = True
  Nt = 9000   # number of timesteps
  
if model == 'narrow':
  animationfile = 'latticenarrow.mp4'  
  
  narrowbasey = 5
  narrowx = [100,100,300]
  narrowheight = [0,45,30]
  narrowwidth = [50,40,90]

  obstacleheight = 30
  obstaclewidth = 50
  tau = 0.6    # collision timescale
  show = 'vorticity'
  cmin = -0.015
  cmax = 0.015
  show = 'velocity'
  cmin = 0
  cmax = 50
  #show = 'density'
  #cmin = -0.01
  #cmax = 0.01  
  plotstream = False
  plotarrow = True
  Nt = 8000   # number of timesteps

if model == 'tube':
  animationfile = 'latticetube_1.mp4'  
  
  tubex = 200
  tuberadius = 75 

  tau = 0.6    # collision timescale
  show = 'vorticity'
  cmin = -0.015
  cmax = 0.015
  show = 'velocity'
  cmin = 0
  cmax = 50
  #show = 'density'
  #cmin = -0.01
  #cmax = 0.01  
  plotstream = False
  plotarrow = True
  Nt = 9000   # number of timesteps

	
# Lattice speeds / weights
NL = 9
idxs = np.arange(NL)
cxs = np.array([0, 0, 1, 1, 1, 0,-1,-1,-1])
cys = np.array([0, 1, 1, 0,-1,-1,-1, 0, 1])
weights = np.array([4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36]) # sums to 1
	
# Initial Conditions
F = np.ones((Ny,Nx,NL)) #* rho0 / NL
np.random.seed(42)
F += 0.01*np.random.randn(Ny,Nx,NL)
X, Y = np.meshgrid(range(Nx), range(Ny))
#F[:,:,3] += 2 #* (1+0.2*np.cos(2*np.pi*X/Nx*4))
F[:,:,3] += 1 #* (1+0.2*np.cos(2*np.pi*X/Nx*4))

rho = np.sum(F,2)
for i in idxs:
	F[:,:,i] *= rho0 / rho

ux  = np.sum(F*cxs,2) / rho
uy  = np.sum(F*cys,2) / rho


def rotate_matrix(matrix, angle):
    # Convert angle to radians
    angle_rad = math.radians(angle)

    # Get matrix dimensions
    rows, cols = len(matrix), len(matrix[0])

    # Get center point
    center_x, center_y = cols / 2, rows / 2

    # Create empty rotated matrix
    rotated_matrix = np.zeros((rows, cols), dtype=bool)

    # Iterate over each pixel in the original matrix
    for y in range(rows):
        for x in range(cols):
            # Translate coordinates to center
            x_trans, y_trans = x - center_x, y - center_y

            # Apply rotation
            x_rotated = x_trans * math.cos(angle_rad) - y_trans * math.sin(angle_rad)
            y_rotated = x_trans * math.sin(angle_rad) + y_trans * math.cos(angle_rad)

            # Translate coordinates back to original position
            x_rotated += center_x
            y_rotated += center_y

            # Interpolate value from original matrix to rotated matrix
            if 0 <= x_rotated < cols - 1 and 0 <= y_rotated < rows - 1:
                x0, y0 = int(x_rotated), int(y_rotated)
                x1, y1 = x0 + 1, y0 + 1

                # Bilinear interpolation
                dx, dy = x_rotated - x0, y_rotated - y0
                top_left = matrix[y0, x0] * (1 - dx) * (1 - dy)
                top_right = matrix[y0, x1] * dx * (1 - dy)
                bottom_left = matrix[y1, x0] * (1 - dx) * dy
                bottom_right = matrix[y1, x1] * dx * dy

                rotated_matrix[y, x] = top_left + top_right + bottom_left + bottom_right
                if rotated_matrix[y, x] >= 0.5: rotated_matrix[y, x] = True
                if rotated_matrix[y, x] < 0.5: rotated_matrix[y, x] = False
                
    return rotated_matrix


#@njit   
def circlecond(X,Y,obstacleradius,obstaclex,obstacley):
    return (X - obstaclex)**2 + (Y - obstacley)**2 < (obstacleradius)**2
def ellipsecond(X,Y,obstacleradiusa,obstacleradiusb,obstaclex,obstacley):
    ell = (X-Nx/2)**2/(obstacleradiusa)**2 + (Y - Ny/2)**2/(obstacleradiusb)**2 < 1
    ell2 = rotate_matrix(ell, stallangle)
    ell2 = np.roll(ell2,-int(Nx/2-circlex[0]),axis=1)
    return ell2
	
# obstacle boundary
X, Y = np.meshgrid(range(Nx), range(Ny))
if model == 'cylinder':
  obstacle = ( circlecond(X,Y,obstacleradius,circlex[0],circley[0]) )
  for i in range(1,len(circlex)):
      obstacle = obstacle | ( circlecond(X,Y,obstacleradius,circlex[i],circley[i]) )
if model == 'square':
  obstacle = ((X >= obstaclex-obstaclewidth/2) & (X <= obstaclex+obstaclewidth/2)) & ((Y <= Ny-(Ny-obstacleheight)/2) & (Y >= (Ny-obstacleheight)/2))
if model == 'narrow':
  obstacle = ((Y > Ny-narrowbasey - narrowheight[0]*np.exp(-(X-narrowx[0])**2/narrowwidth[0]**2)) | 
              (Y < narrowbasey + narrowheight[0]*np.exp(-(X-narrowx[0])**2/narrowwidth[0]**2)))  
  for i in range(1,len(narrowx)):
      obstacle = obstacle | ((Y > Ny - narrowbasey - narrowheight[i]*np.exp(-(X-narrowx[i])**2/narrowwidth[i]**2)) | 
                             (Y < narrowbasey + narrowheight[i]*np.exp(-(X-narrowx[i])**2/narrowwidth[i]**2)))  
if model == 'tube':
  obstacle = (((Y > tuberadius/2+Ny/2) | (Y < -tuberadius/2+Ny/2)) & 
               (X < tubex))  
if model == 'ellipse':
  obstacle = ( ellipsecond(X,Y,obstacleradiusa,obstacleradiusb,circlex[0],circley[0]) )
	
# Prep figure
fig = plt.figure(figsize=(12,6))


vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) - (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1))
vorticity[obstacle] = np.nan
vorticity = np.ma.array(vorticity, mask=obstacle)
im1 = plt.imshow(vorticity, cmap='seismic',vmin = cmin, vmax = cmax)
im2 = plt.imshow(~obstacle, cmap='gray', alpha=0.3)


if plotarrow:
  subsample_factor = 10
  Xsub = X[::subsample_factor, ::subsample_factor]
  Ysub = Y[::subsample_factor, ::subsample_factor]
  uxsub = ux[::subsample_factor, ::subsample_factor]
  uysub = uy[::subsample_factor, ::subsample_factor]
  arrow = plt.quiver(Xsub,Ysub,uxsub,uysub)

if plotstream:
  streams = plt.streamplot(X,Y,ux,uy,density=1,linewidth=0.5,color='b',cmap='jet',arrowsize=0)


ax = plt.gca()
plt.clim(-.1, .1) 
ax.invert_yaxis()
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)	
ax.set_aspect('equal')
if show=='vorticity':	
  cbar = fig.colorbar(im1,ax=ax,label='vorticity')
if show=='density':	
  cbar = fig.colorbar(im1,ax=ax,label='$\Delta$ density')
if show=='velocity':	
  cbar = fig.colorbar(im1,ax=ax,label='velocity')

if model == 'cylinder':
  theta = np.linspace(0,2*np.pi,180)
  for i in range(len(circlex)):
    x = circlex[i]+obstacleradius*np.cos(theta)
    y = circley[i]+obstacleradius*np.sin(theta)
    plt.plot(x,y,color='black',lw=2)

if model == 'ellipse':
  theta = np.linspace(0,2*np.pi,180)
  for i in range(len(circlex)):
    x = obstacleradiusa*np.cos(theta)
    y = obstacleradiusb*np.sin(theta)
    x2 = x*np.cos(stallangle*np.pi/180)+y*np.sin(stallangle*np.pi/180)
    y2 = -x*np.sin(stallangle*np.pi/180)+y*np.cos(stallangle*np.pi/180)
    x2 = x2 + circlex[i]
    y2 = y2 + circley[i]
    plt.plot(x2,y2,color='black',lw=2)

if model == 'square':
  xl = obstaclex-obstaclewidth/2
  xr = obstaclex+obstaclewidth/2
  yo = Ny-(Ny-obstacleheight)/2
  yu = (Ny-obstacleheight)/2
  plt.plot([xl,xl,xr,xr,xl],[yu,yo,yo,yu,yu],color='black',lw=2)

if model == 'narrow':
  x = np.linspace(0,Nx-1,500)
  y1 = Ny-narrowbasey
  y2 = narrowbasey
  for i in range(len(narrowx)):
    y1 = y1 - narrowheight[i]*np.exp(-(x-narrowx[i])**2/narrowwidth[i]**2) 
    y2 = y2 + narrowheight[i]*np.exp(-(x-narrowx[i])**2/narrowwidth[i]**2) 
  plt.plot(x,y1,color='black',lw=2)
  plt.plot(x,y2,color='black',lw=2)
  
if model == 'tube':
  plt.plot([0,tubex,tubex],[Ny/2+tuberadius/2,Ny/2+tuberadius/2,Ny-1],color='black',lw=2)
  plt.plot([0,tubex,tubex],[Ny/2-tuberadius/2,Ny/2-tuberadius/2,0],color='black',lw=2)
  
  
# info box
infobox = ''
infobox += 'collisions time scale tau: ' + "{:.1f}".format(tau) + ''    
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8) 
ax.text(0.02,0.95,infobox, fontsize=10,bbox=props,verticalalignment='top',transform=ax.transAxes)
  
	


#@njit    
#def rollF():
#  global F  
#  for i, cx, cy in zip(idxs, cxs, cys):
#  			F[:,:,i] = np.roll(F[:,:,i], cx, axis=1)
#  			F[:,:,i] = np.roll(F[:,:,i], cy, axis=0)
 
#@njit
#def collisionF():
#  global F  
#  Feq = np.zeros(F.shape)
#  for i, cx, cy, w in zip(idxs, cxs, cys, weights):
#			Feq[:,:,i] = rho * w * ( 1 + 3*(cx*ux+cy*uy)  + 9*(cx*ux+cy*uy)**2/2 - 3*(ux**2+uy**2)/2 )
#		
#  F += -(1.0/tau) * (F - Feq)   

# Simulation Main Loop    
def animate(it):
  global F   
  global arrow, streams     
 		  
  # Drift
  #rollF() 
  for i, cx, cy in zip(idxs, cxs, cys):
  			F[:,:,i] = np.roll(F[:,:,i], cx, axis=1)
  			F[:,:,i] = np.roll(F[:,:,i], cy, axis=0)
		
		
  # Set reflective boundaries
  bndryF = F[obstacle,:]
  bndryF = bndryF[:,[0,5,6,7,8,1,2,3,4]]
	
		
  # Calculate fluid variables
  rho = np.sum(F,2)
  ux  = np.sum(F*cxs,2) / rho
  uy  = np.sum(F*cys,2) / rho
		
		
  # Apply Collision
  # collisionF()
  Feq = np.zeros(F.shape)
  for i, cx, cy, w in zip(idxs, cxs, cys, weights):
  			Feq[:,:,i] = rho * w * ( 1 + 3*(cx*ux+cy*uy)  + 9*(cx*ux+cy*uy)**2/2 - 3*(ux**2+uy**2)/2 )
  		
  F += -(1.0/tau) * (F - Feq)
		
  # Apply boundary 
  F[obstacle,:] = bndryF
	
		
  # plot  color 1/2 particles blue, other half red
  if (it % 1 == 0):
    ux[obstacle] = 0
    uy[obstacle] = 0
    if show=='vorticity':
      vorticity = (np.roll(ux, -1, axis=0) - np.roll(ux, 1, axis=0)) - (np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1))
      vorticity[obstacle] = np.nan
      vorticity = np.ma.array(vorticity, mask=obstacle)
      im1.set_array(vorticity)
    if show=='density':
      density = rho-rho0
      density[obstacle] = np.nan
      density = np.ma.array(100*density, mask=obstacle)
      im1.set_array(density)
    if show=='velocity':
      velocity = np.sqrt(ux**2+uy**2)
      velocity[obstacle] = np.nan
      velocity = np.ma.array(200*velocity, mask=obstacle)
      im1.set_array(velocity)
    
    im2.set_array(~obstacle)
    
    if plotarrow:
      subsample_factor = 10
      Xsub = X[::subsample_factor, ::subsample_factor]
      Ysub = Y[::subsample_factor, ::subsample_factor]
      uxsub = ux[::subsample_factor, ::subsample_factor]
      uysub = uy[::subsample_factor, ::subsample_factor]
      arrow.remove()
      arrow = plt.quiver(Xsub,Ysub,uxsub,uysub,scale=10,alpha=0.5,width=0.002)

    if plotstream:
      streams.lines.remove()
      streams = plt.streamplot(X,Y,ux,uy,density=1,linewidth=0.5,color='b',cmap='jet',arrowsize=0)

    print(it)
  		


anim = animation.FuncAnimation(fig,animate,interval=1,frames=Nt)
anim.save(animationfile,fps=25,dpi=300)


