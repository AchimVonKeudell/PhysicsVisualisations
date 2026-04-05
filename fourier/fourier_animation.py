#
# Fourier Reihen
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation  as animation
import scipy.constants as const


# --------------------------------------------------------------------------
# Choose Model 
# --------------------------------------------------------------------------
model = "Sawtooth"
#model = "Rectangular"
#model = "Triangular"


# --------------------------------------------------------------------------
# Timing Simulation
# --------------------------------------------------------------------------
t = 0
dt = 2 # time step
steps = 720 # simulation steps

# --------------------------------------------------------------------------
# Physical parameter Simulation
# --------------------------------------------------------------------------
  
if model == "Sawtooth":
  animationname = 'Sawtooth.gif'
  name = 'Fourier Series Sawtooth'
  projection = 'x'
  omega = [1,2,3,4,5,6,7,8,9,10]
  radius = [2,1,0.66667,0.5,0.4,0.3333,0.2857,0.25,0.2222,0.2]  
  phase = [90,-90,90,-90,90,-90,90,-90,90,-90]
  modelx = [0,180,180,540,540,900,900,1260,1260,1440]
  modely = [0,-np.pi,np.pi,-np.pi,np.pi,-np.pi,np.pi,-np.pi,np.pi,0]

if model == "Rectangular":
  animationname = 'Rectangular.gif'
  name = 'Fourier Series Rectangular'
  projection = 'y'
  omega = [0,1,2,3,4,5,6,7,8,9]
  radius = [0.5,0.3182,0,0.1061,0,0.0637,0,0.0455,0,0.0354]  
  phase = [0,0,0,0,0,0,0,0,0,0]  
  modelx = [0,180,180,360,360,540,540,720,720,900,900,1080,1080,1260,1260,1440]
  ampl = np.pi/12
  modely = [ampl,ampl,-ampl,-ampl,ampl,ampl,-ampl,-ampl,
            ampl,ampl,-ampl,-ampl,ampl,ampl,-ampl,-ampl]

if model == "Triangular":
  animationname = 'Triangular.gif'
  name = 'Fourier Series Triangular'
  projection = 'x'
  omega = [1,3,5,7,9,11,13,15,17,19]
  phase = [0,0,0,0,0,0,0,0,0,0]
  radius = [0.8106,0.0901,0.0325,0.0165,0.01,0.0067,0.0048,0.0036,0.0028,0.0022]
  modelx = [0,180,360,540,720,900,1080,1260,1440]
  modely = [1,-1,1,-1,1,-1,1,-1,1]
    
  

# --------------------------------------------------------------------------
# Setup Plot
# --------------------------------------------------------------------------
fig, axp = plt.subplots(1,2,figsize=(10,4.5))

axp[1].legend(fontsize=6,loc=1)


# info box
infobox = name+'\n'
infobox += '$n_i$: '
for i in range(len(omega)):
  infobox += "{:.0f}".format(omega[i])+', '
infobox += '\n'
infobox += 'r$_i$: '
for i in range(len(radius)):
  infobox += "{:.4f}".format(radius[i])+', '
infobox += '\n'
infobox += '$\phi_i$: '
for i in range(len(phase)):
  infobox += "{:.0f}".format(phase[i])+'°, '

props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 

fftx = []
ffty1 = []
ffty2 = []
   
   
# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
#    global t,v,x,a
#    global xt,yt,zt
    global t
   
    print("{:.3f}".format(t),' of ',"{:.3f}".format(steps*dt))
    t += dt
    position = [0,0]
    
    axp[0].clear()
    axp[0].set(xlabel="x (AU)",ylabel="y (AU)")    
    axp[0].set(xlim=(-np.sum(radius)*1.1,np.sum(radius)*1.1),ylim=(-np.sum(radius)*1.1,np.sum(radius)*1.1))
    axp[0].axis('off')
    axp[0].text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=axp[0].transAxes)

    # plot circles 
    for i in range(len(omega)):          
      linex = [] 
      liney = []
      for j in range(360):
        linex.append(position[0]+radius[i]*np.cos(j/360*2*np.pi)) 
        liney.append(position[1]+radius[i]*np.sin(j/360*2*np.pi))        
      # trajectory lines
      if projection == 'x':      
        axp[0].plot(liney,linex,color='m',lw=0.5)
      else:
        axp[0].plot(linex,liney,color='m',lw=0.5)
      position[0] = position[0] + radius[i]*np.cos(phase[i]/360*2*np.pi+t/360*omega[i]*2*np.pi)
      position[1] = position[1] + radius[i]*np.sin(phase[i]/360*2*np.pi+t/360*omega[i]*2*np.pi)
    
    # plot marker at the end
    if projection == 'x':  
      rx = [0,position[1]]   
      ry = [0,position[0]]
    else:
      rx = [0,position[0]]   
      ry = [0,position[1]]

    axp[0].plot(rx,ry,color='orange',lw=0.5)
    axp[0].plot(rx[1],ry[1],color='orange',marker='o')
    axp[0].plot([rx[1],rx[1]],[0,ry[1]],color='red',lw=0.5,ls='dashed')

    
    # plot axis    
    axp[0].plot([-np.sum(radius)*1.1,np.sum(radius)*1.1],[0,0],color='black',lw=0.5)
    axp[0].plot([0,0],[-np.sum(radius)*1.1,np.sum(radius)*1.1],color='black',lw=0.5)

    # plot resulting 
    fftx.append(t)
    ffty1.append(position[0])
    ffty2.append(position[1])
    
    axp[1].clear()
    axp[1].set(xlabel="deg (°)",ylabel="amplitude")
    axp[1].set(xlim=(0,steps*dt),ylim=(-np.sum(radius)*1.1,np.sum(radius)*1.1))
    if projection == 'x':
      axp[1].plot(fftx,ffty1,color='r',lw=0.5)
    else:
      axp[1].plot(fftx,ffty2,color='r',lw=0.5)
    axp[1].plot([0,steps*dt],[0,0],color='black',lw=0.5)
    axp[1].plot(modelx,modely,color='black',lw=0.5,ls='dashed') 
    
    # Update the title of the plot
    fig.suptitle('deg: ' + "{:.0f}".format(t)+'°')
  
         
       
# ----------------------------------------
# Create the complete time sequence
# teh total time span is tend/dtplot
# returns the rate in
# ---------------------------------
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)
