#
# Resonance
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.animation as animation
from matplotlib.collections import LineCollection

       
g = 9.8
damping = 0  
dtplot = 0.05
tend = 200
steps = int(tend/dtplot)

# time step has to be set rather small 
# so that numerical error do not accumulate 
# too much, otherwise the phase space 
# will not be covered fully. This is adjusted by
# sublooping = timesteps in between plotting
sublooping = 2000
dt = dtplot/sublooping 
quiverscale = 1

animationfile = 'pendulumseries.gif'
 
# setup plot
xscale = 3

fig, ax = plt.subplots(1,1,figsize=(6,6))
ax.set(xlabel='x (m)',ylabel='y (m)')
ax.set_ylim(-5.5, 0.1)
ax.set_xlim(-xscale, xscale)
ax.set_axis_off()
ax.set_aspect('equal') 
ax.plot([0,0],[-5.5,0.1],lw=0.5,linestyle='dashed')

numberpendulum = 11
l = [2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5]
# Choose a colormap
cmap = plt.get_cmap('viridis')

# Generate evenly spaced values between 0 and 1
colors = [cmap(i) for i in np.linspace(0, 1, numberpendulum)]

pendulum = []
pendulummarker = []
length = []
A1 = []
B1 = []
omega = []
theta = []


for i in range(numberpendulum):

  pendulum.append(ax.plot(0,0,lw=2,color=colors[i])[0])
  pendulummarker.append(ax.plot(0,0, marker='o',markersize=6,lw=0,color=colors[i])[0])
  length.append(l[i])

  A1.append(0.1)
  B1.append(damping)

  omega.append(0)
  theta.append(np.pi*0.2)


t = 0 

# info box
infobox = ''
infobox += 'pendulua: '+"{:.2f}".format(numberpendulum)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.98,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)

     
def animate(k):
    
    global theta,B1,omega
    
    global t,steps

    print('time: '+ "{:.2f}".format(t)+' s of '+"{:.0f}".format(steps*dtplot)+' s ')
    
    t += dt*sublooping
    
    for k in range(numberpendulum):
      for i in range(sublooping):
         domega = -B1[k]*omega[k] - g/length[k]*np.sin(theta[k])
         omega[k] += domega*dt
         theta[k] += omega[k]*dt

    
      pendulum[k].set_xdata([0,-length[k]*np.sin(theta[k])])
      pendulum[k].set_ydata([0,-length[k]*np.cos(theta[k])])    
      pendulummarker[k].set_xdata([-length[k]*np.sin(theta[k])])
      pendulummarker[k].set_ydata([-length[k]*np.cos(theta[k])])   
  
    fig.suptitle('Time: ' + "{:.2f}".format(t)+' s')
        
        
anim = animation.FuncAnimation(fig,animate,interval=0,frames=steps)
anim.save(animationfile,fps=25,dpi=300) 
