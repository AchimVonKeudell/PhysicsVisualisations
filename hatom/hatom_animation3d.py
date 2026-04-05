#
# H atom 3D
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#
# original code by Sebastian Mag | August 2023
# https://github.com/ssebastianmag/hydrogen-wavefunctions

from scipy.constants import physical_constants
import matplotlib.pyplot as plt
import scipy.special as sp
import numpy as np
import matplotlib.animation as animation

model = 'test'
model = 'nupdown'
model = 'mupdown'
model = 'nupmlscan'

modelcmap = 'hot'
modelnb = 5000
modelsize = 0.5
modelalpha = 0.5

showanimation = True

if model == 'test':
  animationfile = 'test3d.gif'
  n = [7]
  l = [4]
  m = [2]
  NBScan = 50
  grid_extent = 3000
if model == 'nupdown':
  animationfile = 'nupdown3d.gif'
  n = [1,2,3,4,5,6,7,7,7,7,7,7,7]
  l = [0,1,2,3,4,5,6,5,4,3,2,1,0]
  m = [0,0,0,0,0,0,0,0,0,0,0,0,0]
  NBScan = 50
  grid_extent = 3000
if model == 'mupdown':
  animationfile = 'mupdown3d.gif'
  n = [7,7,7,7,7,7,7,7,7]
  l = [4,4,4,4,4,4,4,4,4]
  m = [0,1,2,3,4,3,2,1,0]
  NBScan = 50
  grid_extent = 3000
if model == 'nupmlscan':
  animationfile = 'nuplmscan3d.gif'
  n = [1,2,3,4,5,5,5,5,5,5,5,5,5]
  l = [0,1,2,3,4,4,4,4,4,3,2,1,0]
  m = [0,0,0,0,0,1,2,3,4,3,2,1,0]
  NBScan = 50
  grid_extent = 3000 
NBfunctions = len(n)


def radial_function(n, l, r, a0_scale_factor):
    """ Compute the normalized radial part of the wavefunction using
    Laguerre polynomials and an exponential decay factor.

    Args:
        n (int): principal quantum number
        l (int): azimuthal quantum number
        r (numpy.ndarray): radial coordinate
        a0 (float): scaled Bohr radius
    Returns:
        numpy.ndarray: wavefunction radial component
    """
    a0 = a0_scale_factor * physical_constants['Bohr radius'][0] * 1e+12

    laguerre = sp.genlaguerre(n - l - 1, 2 * l + 1)
    p = 2 * r / (n * a0)

    constant_factor = np.sqrt(
        ((2 / n * a0) ** 3 * (sp.factorial(n - l - 1))) /
        (2 * n * (sp.factorial(n + l)))
    )
    return constant_factor * np.exp(-p / 2) * (p ** l) * laguerre(p)


def angular_function_theta(m, l, theta):
    """ Compute the normalized angular part of the wavefunction using
    Legendre polynomials and a phase-shifting exponential factor.

    Args:
        m (int): magnetic quantum number
        l (int): azimuthal quantum number
        theta (numpy.ndarray): polar angle
        phi (int): azimuthal angle
    Returns:
        numpy.ndarray: wavefunction angular component
    """

    legendre = sp.lpmv(m, l, np.cos(theta))

    constant_factor = ((-1) ** m) * np.sqrt(
        ((2 * l + 1) * sp.factorial(l - np.abs(m))) /
        (4 * np.pi * sp.factorial(l + np.abs(m)))
    )
    return constant_factor * legendre

def angular_function_phi(m, l, phi):
    """ Compute the normalized angular part of the wavefunction using
    Legendre polynomials and a phase-shifting exponential factor.

    Args:
        m (int): magnetic quantum number
        l (int): azimuthal quantum number
        theta (numpy.ndarray): polar angle
        phi (int): azimuthal angle
    Returns:
        numpy.ndarray: wavefunction angular component
    """

    return np.real(np.exp(1.j * m * phi))



radius = np.linspace(0,grid_extent,grid_extent)
theta = np.linspace(0,2*np.pi,100)
phi = np.linspace(0,np.pi,100)

wf_rseries = []
wf_thetaseries = []
wf_phiseries = []
wf_cseries = []


for i in range(len(n)):
  print('computing: ',i)  

  wf_radius = []
  wf_theta = []
  wf_phi = []
  wf_c = []

  prob_r = radial_function(n[i], l[i], radius, 1)
  prob_r = np.abs(prob_r)**2
  prob_r = prob_r/np.sum(prob_r)
  
  prob_theta = angular_function_theta(m[i], l[i], theta)
  prob_theta = np.abs(prob_theta)**2
  prob_theta = prob_theta/np.sum(prob_theta)
  
  prob_phi = angular_function_phi(m[i], l[i], phi)
  prob_phi = np.abs(prob_phi)**2
  prob_phi = prob_phi/np.sum(prob_phi)
  

  for k in range(modelnb):
    r = np.random.choice(radius,p=prob_r) 
    th = np.random.choice(theta,p=prob_theta) 
    ph = np.random.choice(phi,p=prob_phi) 
  
    p_r = radial_function(n[i], l[i], r, 1)
    p_theta = angular_function_theta(m[i], l[i], th)
    p_phi = angular_function_phi(m[i], l[i], ph)
  
    wf_c.append(np.abs(p_r*p_theta*p_phi)**2)
    #wf_c.append(1)
    wf_radius.append(r)
    wf_theta.append(th)
    wf_phi.append(ph)
    
  wf_rseries.append(wf_radius)
  wf_thetaseries.append(wf_theta)
  wf_phiseries.append(wf_phi)  
  wf_cseries.append(wf_c)
  

#fig, ax = plt.subplots()

fig = plt.figure(figsize=(8, 8.5))
ax = fig.add_subplot(111, projection='3d')
plt.tight_layout()

x = wf_rseries[0]*np.sin(wf_thetaseries[0])*np.cos(wf_phiseries[0])
y = wf_rseries[0]*np.sin(wf_thetaseries[0])*np.sin(wf_phiseries[0])
z = wf_rseries[0]*np.cos(wf_thetaseries[0])
c = wf_cseries[0]

ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
ax.set_axis_off()
ax.plot([0,0],[0,0],[-grid_extent,grid_extent],color='red')
ax.plot([0,0],[-grid_extent,grid_extent],[0,0],color='red')
ax.plot([-grid_extent,grid_extent],[0,0],[0,0],color='red')
  
im = ax.scatter3D(x, y, z, c=c, s=modelsize, cmap=modelcmap,alpha=modelalpha)
cbar = fig.colorbar(im, ax = ax, shrink = 0.4, aspect = 20)
cbar.set_ticks([])
 


def animate(k):
 # Compute and visualize the wavefunction probability density
    
    #global xb, yb
    global grid_extent
    
    i = int(k/NBScan)
    print(k,': ',n[i],',',l[i],',',m[i])
    NBalpha = k-i*NBScan
    
    ax.clear()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    #ax.set_xlim(-grid_extent*a0,grid_extent*a0)
    #ax.set_ylim(-grid_extent*a0,grid_extent*a0)
    ax.set_axis_off()
    
    fig.suptitle('(n,l,m): ('+"{:.0f}".format(n[i])+','
                               +"{:.0f}".format(l[i])+','
                               +"{:.0f}".format(m[i])+')',fontsize=20)

    
    x = wf_rseries[i]*np.sin(wf_thetaseries[i])*np.cos(wf_phiseries[i])
    y = wf_rseries[i]*np.sin(wf_thetaseries[i])*np.sin(wf_phiseries[i])
    z = wf_rseries[i]*np.cos(wf_thetaseries[i])
    c = wf_cseries[i]
    
    
    if NBalpha < NBScan/2:
      ax.scatter3D(x, y, z, c=c, s=modelsize, cmap=modelcmap,alpha=NBalpha/(NBScan/2)*modelalpha)
    else:
      ax.scatter3D(x, y, z, c=c, s=modelsize, cmap=modelcmap,alpha=(1-(NBalpha-NBScan/2)/(NBScan/2))*modelalpha)
    
    ax.plot([0,0],[0,0],[-grid_extent,grid_extent],color='red')
    ax.plot([0,0],[-grid_extent,grid_extent],[0,0],color='red')
    ax.plot([-grid_extent,grid_extent],[0,0],[0,0],color='red')
    
    ax.view_init(elev=60*np.sin(k/(NBScan*NBfunctions*0.5)*np.pi)**2, azim=45+k)
    
    
if showanimation:
  anim = animation.FuncAnimation(fig,animate,interval=1,frames=NBScan*NBfunctions)
  anim.save(animationfile,fps=25,dpi=300)
