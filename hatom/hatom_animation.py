#
# H atom
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
#model = 'nupdown'
#model = 'mupdown'
model = 'nupmlscan'

if model == 'test':
  animationfile = 'test.gif'
  n = [1]
  l = [0]
  m = [0]
  NBScan = 2
  grid_extent = 3000
if model == 'nupdown':
  animationfile = 'nupdown.gif'
  n = [1,2,3,4,5,6,7,7,7,7,7,7,7]
  l = [0,1,2,3,4,5,6,5,4,3,2,1,0]
  m = [0,0,0,0,0,0,0,0,0,0,0,0,0]
  NBScan = 50
  grid_extent = 3000
if model == 'mupdown':
  animationfile = 'mupdown.mp4'
  n = [7,7,7,7,7,7,7,7,7]
  l = [4,4,4,4,4,4,4,4,4]
  m = [0,1,2,3,4,3,2,1,0]
  NBScan = 50
  grid_extent = 7000
if model == 'nupmlscan':
  animationfile = 'nuplmscan.mp4'
  n = [1,2,3,4,5,5,5,5,5,5,5,5,5]
  l = [0,1,2,3,4,4,4,4,4,3,2,1,0]
  m = [0,0,0,0,0,1,2,3,4,3,2,1,0]
  NBScan = 50
  grid_extent = 3000 



def radial_function(n, l, r, a0):
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

    laguerre = sp.genlaguerre(n - l - 1, 2 * l + 1)
    p = 2 * r / (n * a0)

    constant_factor = np.sqrt(
        ((2 / n * a0) ** 3 * (sp.factorial(n - l - 1))) /
        (2 * n * (sp.factorial(n + l)))
    )
    return constant_factor * np.exp(-p / 2) * (p ** l) * laguerre(p)


def angular_function(m, l, theta, phi):
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
    return constant_factor * legendre * np.real(np.exp(1.j * m * phi))


def compute_wavefunction(n, l, m, a0_scale_factor):
    """ Compute the normalized wavefunction as a product
    of its radial and angular components.

    Args:
        n (int): principal quantum number
        l (int): azimuthal quantum number
        m (int): magnetic quantum number
        a0_scale_factor (float): Bohr radius scale factor
    Returns:
        numpy.ndarray: wavefunction
    """

    # Scale Bohr radius for effective visualization
    a0 = a0_scale_factor * physical_constants['Bohr radius'][0] * 1e+12

    # z-x plane grid to represent electron spatial distribution
    grid_resolution = 680
    z = x = np.linspace(-grid_extent, grid_extent, grid_resolution)
    z, x = np.meshgrid(z, x)

    # Use epsilon to avoid division by zero during angle calculations
    eps = np.finfo(float).eps

    # Ψnlm(r,θ,φ) = Rnl(r).Ylm(θ,φ)
    psi = radial_function(
        n, l, np.sqrt((x ** 2 + z ** 2)), a0
    ) * angular_function(
        m, l, np.arctan(x / (z + eps)), 0
    )
    return psi

def compute_bohrradius(a0_scale_factor):
    """ Compute the normalized wavefunction as a product
    of its radial and angular components.

    Args:
        n (int): principal quantum number
        l (int): azimuthal quantum number
        m (int): magnetic quantum number
        a0_scale_factor (float): Bohr radius scale factor
    Returns:
        numpy.ndarray: wavefunction
    """

    # Scale Bohr radius for effective visualization
    a0 = a0_scale_factor * physical_constants['Bohr radius'][0] * 1e+12

    theta = np.linspace(0,np.pi*2,360)
    print(grid_extent, a0)
    x = a0*np.cos(theta)
    y = a0*np.sin(theta)
    
    return x, y


def compute_probability_density(psi):
    """ Compute the probability density of a given wavefunction.
    Args:
        psi (numpy.ndarray): wavefunction
    Returns:
        numpy.ndarray: wavefunction probability density
    """
    return np.abs(psi) ** 2


prob_density = []

for i in range(len(n)):
  psi = compute_wavefunction(n[i], l[i], m[i], 1)
  prob_density.append(compute_probability_density(psi))

NBfunctions = len(n)

fig, ax = plt.subplots(figsize=(16, 16.5))
im = ax.imshow(np.sqrt(prob_density[0]).T, cmap='Reds')
cbar = plt.colorbar(im, fraction=0.046, pad=0.03)
cbar.set_ticks([])
#a0 = physical_constants['Bohr radius'][0] * 1e+12
#xb, yb = compute_bohrradius(1)


def animate(k):
 # Compute and visualize the wavefunction probability density
    
    #global xb, yb
    
    i = int(k/NBScan)
    print(k,': ',n[i],',',l[i],',',m[i])
    NBalpha = k-i*NBScan
    
    ax.clear()
    ax.set_xticks([])
    ax.set_yticks([])
    #ax.set_xlim(-grid_extent*a0,grid_extent*a0)
    #ax.set_ylim(-grid_extent*a0,grid_extent*a0)
    ax.text(30,50,'(n,l,m): ('+"{:.0f}".format(n[i])+','
                               +"{:.0f}".format(l[i])+','
                               +"{:.0f}".format(m[i])+')',fontsize=40)
    if NBalpha < NBScan/2:
      ax.imshow(np.sqrt(prob_density[i]).T, cmap='Reds',alpha=NBalpha/(NBScan/2))
    else:
      ax.imshow(np.sqrt(prob_density[i]).T, cmap='Reds',alpha=1-(NBalpha-NBScan/2)/(NBScan/2))
    
    ax.text(30,-30,'H atom, electron probabillty density',fontsize=30)
    #ax.plot(xb,yb)
    
anim = animation.FuncAnimation(fig,animate,interval=1,frames=NBScan*NBfunctions)
anim.save(animationfile,fps=25,dpi=300)
