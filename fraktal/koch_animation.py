#
# Fractrals Koch Curve
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2026
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation  as animation

model = 'zoom'
#model = 'statisch'

if model =='zoom':
  animationname = 'koch_zoom.gif'
if model =='statisch':
  animationname = 'koch_dynamisch.gif'

faktor = -np.sqrt(3) / 2
x = [-80,0,80,-80]
y = [-80*(-faktor)+25,80*(-faktor)+25,-80*(-faktor)+25,-80*(-faktor)+25]


def koch(x1,x2,y1,y2):

  xv = []
  yv = []
  
  # line 1
  xv.append(x1)
  xv.append((2 * x1 + x2) / 3)
  yv.append(y1)
  yv.append((2 * y1 + y2) / 3)
  
  # line 2
  xv.append((2 * x1 + x2) / 3)
  xv.append((x1 + x2) / 2 + faktor * (y2 - y1) / 3)
  yv.append((2 * y1 + y2) / 3)
  yv.append((y1 + y2) / 2 + faktor * (x1 - x2) / 3)
  
  # line 3
  xv.append((x1 + x2) / 2 + faktor * (y2 - y1) / 3)
  xv.append((x1 + 2 * x2) / 3)
  yv.append((y1 + y2) / 2 + faktor * (x1 - x2) / 3)
  yv.append((y1 + 2 * y2) / 3)
  
  # line 4
  xv.append((x1 + 2 * x2) / 3)
  xv.append(x2)
  yv.append((y1 + 2 * y2) / 3)
  yv.append(y2)
  return xv,yv


# Create Koch Curves:
iterations = 8
kochcurvex = []
kochcurvey = []
kochcurvex.append(x)
kochcurvey.append(y)

def gradient_colors(n, cmap_name="viridis"):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(i) for i in np.linspace(0, 1, n)]

colormap = gradient_colors(iterations, "plasma")

for i in range(iterations):
    xn = []
    yn = []
    for j in range(len(x)-1):
      xv,yv = koch(x[j],x[j+1],y[j],y[j+1])    
      xn += xv
      yn += yv

    x = xn
    y = yn
    kochcurvex.append(x)
    kochcurvey.append(y)
 

fig, ax = plt.subplots(1, 1, figsize=(5, 5),
                        layout="constrained")
  
sweep = 80


# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
   
  global x,y  
  
  kindex = int(k/sweep)  
  print(kindex,' von ',iterations)  

  ax.clear()
  ax.set_aspect('equal')
  ax.axis('off')
  if model == 'zoom':
    move = k/(iterations*sweep)
  if model == 'statisch':
    move = 0
  x0 = 0 - move*65 
  y0 = 0
  xscale = 5+95*(1-move)
  yscale = 5+95*(1-move)
  xmin = x0-xscale
  xmax = x0+xscale
  ymin = y0-yscale
  ymax = y0+yscale 
  ax.set(xlim=(xmin,xmax),ylim=(ymin,ymax))

  # plot current iteration
  ax.plot(kochcurvex[kindex],kochcurvey[kindex],color=colormap[kindex],alpha=(k-kindex*sweep)/(sweep))
  # fade last iteration
  if kindex>0:  
    ax.plot(kochcurvex[kindex-1],kochcurvey[kindex-1],color=colormap[kindex-1],alpha=1-(k-kindex*sweep)/(sweep))
    
  fig.suptitle('Koch Snowflake, Iteration: '+"{:.0f}".format(kindex))

anim = animation.FuncAnimation(fig,animate,interval=1,frames=iterations*sweep)
anim.save(animationname,fps=25,dpi=300)   