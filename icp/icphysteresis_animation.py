# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:03:53 2024

@author: Achim
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap



ne = np.linspace(0,3e11,1000)

nkap = 1e8
nind = 1e11
Irf = 30
Rabs = 30

def pabs(neel):
  return 0.5*Irf*Rabs*(nind*neel/(nind**2+neel**2)+nkap/(nkap+neel))

def ploss(neel):
  return 1e-9*neel+30

def balance(neel):
  return abs(pabs(neel)-ploss(neel))  


cross = np.ones(4)

fig, ax = plt.subplots(1,2,figsize=(8,4))

ax[0].set(xlim=(0,30),ylim=(0,700),xlabel='n$_e$ (x 10$^{10}$ cm$^{-3}$)',ylabel='P (W)')
line1 = ax[0].plot(ne,pabs(ne),label='P$_{abs}$')[0]
line2 = ax[0].plot(ne,ploss(ne),label='P$_{loss}$')[0]
line3 = ax[0].plot(cross,ploss(cross),marker='o',color='red',markersize=5,lw=0)[0]
#ax[0].set_xticklabels([])
#ax[0].set_yticklabels([])
ax[0].legend(loc=1)


linex = []
liney = []
line4 = ax[1].plot(linex,liney)[0]
line5 = ax[1].plot(linex,liney,marker='o',color='red',markersize=5,lw=0)[0]
ax[1].set(xlim=(0,100),ylim=(0,500),xlabel='I$_{rf}$ (A)',ylabel='P (W)')
#ax[1].set_xticklabels([])
#ax[1].set_yticklabels([])

# Create a meshgrid for gradient coloring
# Create a gradient from red to white
gradient = np.linspace(0, 1, 512)#.reshape(128, 1)  # Creates a linear gradient
gradient = np.vstack((gradient, gradient))  # Makes it 2D for displaying as an image

# Display the gradient using imshow
im = ax[1].imshow(gradient, extent=[0.05, 0.07, 0.01, 0.02], origin='lower', cmap='Reds', aspect='auto')

# Define the rectangle parameters: (x, y, width, height)
rectangle = plt.Rectangle((18, 340), 2, 120, edgecolor='lightblue', facecolor='lightblue', linewidth=2)
ax[1].add_patch(rectangle)
rectangle = plt.Rectangle((50, 340), 2, 120, edgecolor='lightblue', facecolor='lightblue', linewidth=2)
ax[1].add_patch(rectangle)
ax[1].text(10,470,'ground')
ax[1].text(42,470,'window')
ax[1].text(28,320,'plasma')



# info box
infobox = ''
infobox += 'n$_{kap}$: ' + "{:.0e}".format(nkap) + ' (cm$^{-3}$)\n' 
infobox += 'n$_{ind}$: ' + "{:.0e}".format(nind) + ' (cm$^{-3}$)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[0].text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax[0].transAxes)



rampNB = 300
Iramp = 80

def animate(k):
    global Irf
    global im
    
    if k<rampNB:
      Irf = k/rampNB*Iramp
    else:
      Irf = Iramp-(k-rampNB)/rampNB*Iramp

    NBStart = 20
    cross = np.ones(NBStart)*1e12
    for i in range(NBStart):
      netest=1e8+i/NBStart*1e11    
      result = fsolve(balance,netest)   
      if ((balance(result[0])<1) and (result[0]>0)): cross[i]=result[0]
    
    #print(cross)
        
    if k<rampNB:
       nex = min(cross)
    else:
       # remove maxima from the plot 
       for i in range(NBStart):
         if cross[i]>9e11: cross[i]=-1e9  
       nex = max(cross)
        
    if nex!=1e12: # and ploss(nex)>0:   
      linex.append(Irf)
      liney.append(ploss(nex))
      line5.set_xdata(Irf)
      line5.set_ydata(ploss(nex))    
    line4.set_xdata(linex)
    line4.set_ydata(liney)
    
    #  print(cross)
    line1.set_xdata(ne/1e10)
    line1.set_ydata(pabs(ne))
    line2.set_xdata(ne/1e10)
    line2.set_ydata(ploss(ne))
    line3.set_xdata(cross/1e10)
    line3.set_ydata(pabs(cross))
    
    lamb = 20*1/(nex/1e11)
    if lamb>30: lamb = 30
    if lamb<0: lamb = 0
    
    im.remove()
    im = ax[1].imshow(gradient, extent=[50-lamb, 50, 350, 450], origin='lower', cmap='Reds', aspect='auto')



    
anim = animation.FuncAnimation(fig,animate,interval=1200,frames=2*rampNB)
anim.save('icphysteresis.mp4',fps=25,dpi=300)    