#
# plasma in a 1d1v PIC code
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#model = 'prop. pot'
#model = 'prop. pot**2'
model = 'prop. sqrt(pot)'
#model = 'prop. exp(-act. energ/pot)'
#model = 'prop. pot*exp(-act. energ/pot)'

scalefactor = 1
Xmax = scalefactor*300
Ymax = Xmax

if model == 'prop. pot':
  modelID = 1  
  activationenergy = 5
  animationname = 'dbm_propot.mp4'
  steps = 3000
if model == 'prop. pot**2':
  modelID = 2  
  activationenergy = 5
  animationname = 'dbm_proppot2.mp4'
  steps = 3000
if model == 'prop. sqrt(pot)':
  modelID = 3  
  activationenergy = 5
  animationname = 'dbm_propsqrtpot.mp4'
  steps = 2000
if model == 'prop. exp(-act. energ/pot)':
  modelID = 4  
  activationenergy = 5
  animationname = 'dbm_propexp.mp4'
  steps = 300
if model == 'prop. pot*exp(-act. energ/pot)':
  modelID = 5  
  activationenergy = 5
  animationname = 'dbm_proppotexp.mp4'
  steps = 200

radiuselectrode = 2

growthsite = np.zeros((Xmax,Ymax))
pattern = np.zeros((Xmax,Ymax))
patternx = []
patterny = []

potential = np.zeros((Xmax,Ymax))
potentialnew = np.zeros((Xmax,Ymax))

radius = Xmax/2.1 # position outer boundary
radiussimulation = radius*0.8 # stop simulation if this is reached

# Set Boundary Condition 
radius = Xmax/2.1 
x0 = int(Xmax/2)
y0 = int(Ymax/2)

# Central Starting point 
pattern[x0,y0]=1
potential[x0,y0]=0

# Inner Electrode 
for i in range(360):
      angle = i/360*(2*np.pi)
      xp = int(np.cos(angle)*radiuselectrode+x0)
      yp =int(np.sin(angle)*radiuselectrode+y0)
      potential[xp,yp] =0
      pattern[xp,yp] =1
      patternx.append(xp)
      patterny.append(yp)

# Outer Boundary 
for i in range(360):
      angle = i/360*(2*3.141);
      xp = int(np.cos(angle)*radius+x0);
      yp = int(np.sin(angle)*radius+y0);
      potential[xp,yp] = scalefactor


# --------------------------------------
#  Solve Poisson Equation
#   Gauss Seidel method
# --------------------------------------}
def SolvePoisson(k):
    poissoniteration = 0
    poissoniterationlimit = 1000
    while poissoniteration<poissoniterationlimit:
        poissoniteration += 1
        # Gauss Seidel step
        for i in range(2,Xmax-1):
            for j in range(2,Ymax-1):
                potentialnew[i,j] = 0.25*(potential[i-1,j]+potential[i+1,j]
                                          +potential[i,j-1]+potential[i,j+1])
         
        # Fix Pattern Condition 
        for i in range(1,Xmax):
           for j in range(1,Ymax):
            if pattern[i,j]==1 : potentialnew[i,j] =0
   
        # Set Outer Boundary Condition  
        x0 = int(Xmax/2)
        y0 = int(Ymax/2)
        for i in range(360):
            angle = i/360*(2*np.pi)
            xp = int(np.cos(angle)*radius+x0)
            yp = int(np.sin(angle)*radius+y0)
            potentialnew[xp,yp] = scalefactor
            potentialnew[xp,yp] = scalefactor
            potentialnew[x0,y0] = 0
  
        # Set Electrode Boundary Condition  
        for i in range(360):
            angle = i/360*(2*np.pi);
            xp = int(np.cos(angle)*radiuselectrode+x0);
            yp = int(np.sin(angle)*radiuselectrode+y0);
            potential[xp,yp] =0
            potential[xp,yp] =0
            pattern[xp,yp] =1
         
        # write array back and check for deviation
        diffmax = 0;
        for i in range(Xmax):
            for j in range(Ymax):
               if abs(potential[i,j]-potentialnew[i,j])>diffmax:
                   diffmax = abs(potential[i,j]-potentialnew[i,j])
               potential[i,j] = potentialnew[i,j]
        if diffmax < 1/Xmax**2*10: return # return when diff < threshold
        if k == poissoniterationlimit: return # return when max iter. reached
     

#----------------------------------
# Define Plot
#----------------------------------
#figsize=(4.5,4.5)
fig = plt.figure(constrained_layout=True)
ax = []
ax.append(fig.add_subplot())
ax[0].set(xlim=(0,Xmax),ylim=(0,Ymax),xlabel='x',ylabel='y')
im = plt.imshow(potential,cmap='YlOrBr',label='potential')
scatter = plt.plot(patternx,patterny,color='r',lw=0,marker='s',markersize=0.5,label='pattern')[0]
fig.colorbar(im,ax=ax[0],label='potential')
ax[0].legend(fontsize=6)

# info box
infobox = ''
infobox += 'model: ' + model + '\n' 
infobox += 'act. energy: ' + "{:.1f}".format(activationenergy) + ''    
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8) 
ax[0].text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax[0].transAxes)

#SolvePoisson(1)
#np.savetxt('potinitial300.dat', np.transpose(potential), fmt='%.3f', header="x1 v1")
potential = np.loadtxt('potinitial300.dat')

# ----------------------------------
# Main loop
# ----------------------------------
def animate(k):

    # Gauss Seidel solution
    SolvePoisson(k)
   
    # Identify Growth Sites
    for i in range(Xmax):
      for j in range(Ymax): growthsite[i,j] = 0
    for i in range(2,Xmax-1):
        for j in range(2,Ymax-1):
           # Check if i,j are within a given circle
           # around the origin
           if pattern[i,j]==0:
                distance = np.sqrt(((i-x0)**2+(j-y0)**2));
                if distance < radiussimulation:
                    if pattern[i+1,j]==1: growthsite[i,j]=1;
                    if pattern[i-1,j]==1: growthsite[i,j]=1;
                    if pattern[i,j+1]==1: growthsite[i,j]=1;
                    if pattern[i,j-1]==1: growthsite[i,j]=1;
  
    # Normalization Growth Site probailities
    GrowthSitesNB = 0
    sum = 0
    for i in range(Xmax):
      for j in range(Ymax):
        if growthsite[i,j]==1:
            GrowthSitesNB += 1
            if modelID == 1:
               sum += potential[i,j]
            elif modelID == 2:   
               sum += (potential[i,j])**2
            elif modelID == 3:   
               sum += np.sqrt(potential[i,j])
            elif modelID == 4:   
               if potential[i,j]!=0: 
                   sum += np.exp(-activationenergy/potential[i,j])
            elif modelID == 5:      
               if potential[i,j]!=0:
                   sum += potential[i,j]*np.exp(-activationenergy/potential[i,j])
            else:
              sum += potential[i,j]
    
    #  Expand Pattern depending on Potential
    for i in range(Xmax):
        for j in range(Ymax):
          weight = 0  
          if growthsite[i,j]==1:
              if modelID == 1:
                weight = potential[i,j]/sum
              elif modelID == 2:  
                weight = (potential[i,j])**2/sum;
              elif modelID == 3:  
                weight = np.sqrt(potential[i,j])/sum;
              elif modelID == 4:  
                if (potential[i,j]!=0) and (sum !=0):
                    weight = np.exp(-activationenergy/potential[i,j])/sum
              elif modelID == 5:  
                if (potential[i,j]!=0) and (sum !=0):
                    weight = potential[i,j]*np.exp(-activationenergy/potential[i,j])/sum
              else:
                if sum!=0:
                    weight = potential[i,j]/sum
                else: weight = 0
              if np.random.rand()<weight:
                pattern[i,j] = 1
                patternx.append(j)
                patterny.append(i)

    print('Step ',k,' of ',steps)
    im.set_array(potential)
    scatter.set_xdata(patternx)
    scatter.set_ydata(patterny)
 
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=180)