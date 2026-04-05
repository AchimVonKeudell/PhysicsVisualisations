#
# Barnsley Fern Fractal
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2026
#

import random
import matplotlib.pyplot as plt
import matplotlib.animation  as animation

mode = 'original'
mode = 'mutant'
#mode = 'mutant2'

x, y = 0.0, 0.0
iterations = 1000
subiterations = 100

if mode == 'original':
  animationname = 'barnsleyfern_original.gif'
  name = 'Orginal Barnsley Fern'
  f1 = [    0,     0,     0, 0.16,  0,    0, 0.01]
  f2 = [ 0.85,  0.04, -0.04, 0.86,  0, 1.60, 0.85]
  f3 = [ 0.20, -0.26,  0.23, 0.22,  0, 1.60, 0.07]
  f4 = [-0.15,  0.28,  0.26, 0.24,  0, 0.44, 0.07]

if mode == 'mutant':
  animationname = 'barnsleyfern_mutant.gif'
  name = 'Orginal Barnsley Mutant 1'
  f1 = [     0,     0,      0, 0.25,      0,  0.4, 0.02]
  f2 = [  0.95, 0.005, -0.005, 0.93, -0.002,  0.5, 0.84]
  f3 = [ 0.035,  -0.2,   0.16, 0.04,  -0.09, 0.02, 0.07]
  f4 = [ -0.04,   0.2,   0.16, 0.04,  0.083, 0.12, 0.07]

if mode == 'mutant2':
  animationname = 'barnsleyfern_mutant2.gif'
  name = 'Orginal Barnsley Mutant 2'
  f1 = [     0,     0,     0, 0.25,  0, -0.14, 0.02]
  f2 = [  0.85,  0.02, -0.02, 0.83,  0,     1, 0.84]
  f3 = [  0.09, -0.28,   0.3, 0.11,  0,   0.6, 0.07]
  f4 = [ -0.09,  0.28,   0.3, 0.09,  0,   0.7, 0.07]


fig, ax = plt.subplots(1, 1, figsize=(8, 8),
                        layout="constrained")

ax.set_aspect('equal')
ax.axis('off')

ax.axis("off")
# info box
infobox = name+'\n'
infobox += 'Coefficients: \n'
infobox +=('f1: ' + "{:.0f}".format(f1[0]) + ', ' + "{:.3f}".format(f1[1]) + ', '
                               + "{:.3f}".format(f1[2]) + ', ' + "{:.3f}".format(f1[3]) + ', '
                               + "{:.3f}".format(f1[4]) + ', ' + "{:.3f}".format(f1[5]) + ', '
                               + "{:.3f}".format(f1[6]) + '\n')
infobox +=('f2: ' + "{:.3f}".format(f2[0]) + ', ' + "{:.3f}".format(f2[1]) + ', '
                               + "{:.3f}".format(f2[2]) + ', ' + "{:.3f}".format(f2[3]) + ', '
                               + "{:.3f}".format(f2[4]) + ', ' + "{:.3f}".format(f2[5]) + ', '
                               + "{:.3f}".format(f2[6]) + '\n')
infobox +=('f3: ' + "{:.3f}".format(f3[0]) + ', ' + "{:.3f}".format(f3[1]) + ', '
                               + "{:.3f}".format(f3[2]) + ', ' + "{:.3f}".format(f3[3]) + ', '
                               + "{:.3f}".format(f3[4]) + ', ' + "{:.3f}".format(f3[5]) + ', '
                               + "{:.3f}".format(f3[6]) + '\n')
infobox +=('f4: ' + "{:.3f}".format(f4[0]) + ', ' + "{:.3f}".format(f4[1]) + ', '
                               + "{:.3f}".format(f4[2]) + ', ' + "{:.3f}".format(f4[3]) + ', '
                               + "{:.3f}".format(f4[4]) + ', ' + "{:.3f}".format(f4[5]) + ', '
                               + "{:.3f}".format(f4[6]))


props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)


def animate(k):

  global x,y
      
  print(k, 'of ',iterations)  
  

  xs, ys = [], []
  for _ in range(subiterations):
        r = random.random()

    #    if r < 0.01:
    #        x, y = 0.0, 0.16 * y
    #    elif r < 0.86:
    #        x, y = 0.85 * x + 0.04 * y, -0.04 * x + 0.85 * y + 1.6
    #    elif r < 0.93:
    #        x, y = 0.20 * x - 0.26 * y, 0.23 * x + 0.22 * y + 1.6
    #    else:
    #        x, y = -0.15 * x + 0.28 * y, 0.26 * x + 0.24 * y + 0.44


        if r < f1[6]:
            x, y = f1[2], f1[3] * y
        elif r < f1[6]+f2[6]:
            x, y = f2[0] * x + f2[1] * y + f2[4], f2[2] * x + f2[3] * y + f2[5]
        elif r < f1[6]+f2[6]+f3[6]:
            x, y = f3[0] * x +f3[1] * y  + f3[4], f3[2] * x + f3[3] * y + f3[5]
        else:
            x, y = f4[0] * x + f4[1] * y + f4[4], f4[2] * x + f4[3] * y + f4[5]

        xs.append(x)
        ys.append(y)

  ax.scatter(xs,ys,s=0.15,color='b')
  fig.suptitle('Barnsley Fern, Random Points: '+"{:.0f}".format(k*1000))

anim = animation.FuncAnimation(fig,animate,interval=1,frames=iterations)
anim.save(animationname,fps=25,dpi=300)   

