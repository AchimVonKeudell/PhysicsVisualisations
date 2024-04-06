#
#
# 
#
# Ising Modell

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numba import jit


framesNB = 2000
microsteps = 50000

T1 = 2.4  
J = 1 
H = 0
animationname = 'ising2T_1.mp4'
T2 = 1  
    
# initialize lattice
size = 200
lattice1 = np.zeros((size,size))
lattice2 = np.zeros((size,size))
for i in range(size):
    for j in range(size):
        if np.random.rand()>0.5:
           lattice1[i,j] = 1
        else:
           lattice1[i,j] = -1 
        if np.random.rand()>0.5:
           lattice2[i,j] = 1
        else:
           lattice2[i,j] = -1 


# setup plot
fig, ax = plt.subplots(1,2,figsize=(8,4.5))
im1 = ax[0].imshow(lattice1,cmap='plasma')
im2 = ax[1].imshow(lattice2,cmap='plasma')
ax[0].set(xlabel="x",ylabel="y")
ax[1].set(xlabel="x",ylabel="y")
fig.suptitle('Ising model, T$_c$ = 2.27')
 

# setup info box
infobox = ''
#infobox += 'J: ' + "{:.2f}".format(J) + ' \n' 
#infobox += 'H: ' + "{:.2f}".format(H) + ' \n' 
infobox += 'T: ' + "{:.1f}".format(T1) + ' (K) ' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=1) 
ax[0].text(0.02,0.98,infobox, fontsize=10,bbox=props,verticalalignment='top',transform=ax[0].transAxes)

# setup info box
infobox = ''
#infobox += 'J: ' + "{:.2f}".format(J) + ' \n' 
#infobox += 'H: ' + "{:.2f}".format(H) + ' \n' 
infobox += 'T: ' + "{:.1f}".format(T2) + ' (K) ' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=1) 
ax[1].text(0.02,0.98,infobox, fontsize=10,bbox=props,verticalalignment='top',transform=ax[1].transAxes)



@jit(nopython = True)
def evolve(lattice,T):
  for i in range(microsteps):
    # spin flip
    x = int(np.random.rand()*size)
    y = int(np.random.rand()*size)
    
    s0 = lattice[x,y]
    
    # calulate current energy    
    xp = x +1
    if xp>size-1: xp = 0
    a1 = s0*lattice[xp,y]
    xm = x -1
    if xm<0: xm = size-1
    a2 = s0*lattice[xm,y]
    yp = y +1
    if yp>size-1: yp = 0
    a3 = s0*lattice[x,yp]
    ym = y -1
    if ym<0: ym = size-1
    a4 = s0*lattice[x,ym]    
    energy1 = -J*(a1+a2+a3+a4)-H*s0
    
    # calulate spin flip
    s0 = (-1)*s0
    xp = x +1
    if xp>size-1: xp = 0
    a1 = s0*lattice[xp,y]
    xm = x -1
    if xm<0: xm = size-1
    a2 = s0*lattice[xm,y]
    yp = y +1
    if yp>size-1: yp = 0
    a3 = s0*lattice[x,yp]
    ym = y -1
    if ym<0: ym = size-1
    a4 = s0*lattice[x,ym]    
    energy2 = -J*(a1+a2+a3+a4)-H*s0
    
    deltae = energy2-energy1
    if deltae <0: # accept flip if energy is more favorable
        lattice[x,y] = lattice[x,y]*(-1)
    else:  # otherwise use activation energy to decide 
        p = np.exp(-deltae/T)
        if np.random.rand()<p:
           lattice[x,y] = lattice[x,y]*(-1)
  
def animate(k):
  
  global lattice1, lattice2
  
  print('frame: ',k,' of ',framesNB)  
  evolve(lattice1,T1)  
  evolve(lattice2,T2)  
  im1.set_array(lattice1)
  im2.set_array(lattice2)
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=framesNB)
anim.save(animationname,fps=25,dpi=300)