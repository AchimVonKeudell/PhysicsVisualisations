#
# Sierpinski Tree
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2026
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation  as animation

model = 'straight'
model = 'curved'
model = 'bending'
model = 'straightrandom'

if model =='straight':
  name = 'Straight Tree'
  generations = 12
  angleleft = 30
  angleright = 30
  branchscale = 0.7
  sweep = 20
  randombranch = 0
  randomangle = 0
  animationname = 'sierpinski_straight_zoom.gif'
if model =='straightrandom':
  name = 'Straight Tree, Random Branches'
  generations = 14
  angleleft = 30
  angleright = 30
  branchscale = 1  
  sweep = 20
  randombranch = 0.3
  randomangle =15
  animationname = 'sierpinski_straight_random_zoom4.gif'
if model =='curved':
  name = 'Curved Tree'  
  generations = 12
  angleleft = 10
  angleright = 30
  branchscale = 0.7
  sweep = 20
  randombranch = 0
  randomangle = 0
  animationname = 'sierpinski_curved_zoom.gif'
if model =='bending':
  name = 'Bending Tree'  
  generations = 12
  angleleft = 5
  angleright = 45
  branchscale = 0.7
  sweep = 20
  randombranch = 0
  randomangle = 0
  animationname = 'sierpinski_bending_zoom.gif'

def generate_sierpinski_tree(
    generations=6,
    initial_length=1.0,
    anglel_deg=30,
    angler_deg=30,
    scale=0.5
):
    """
    Generate Sierpinski tree branches iteratively (no recursion).

    Each branch is represented as:
    (x_start, y_start, length, angle_in_radians)
    """

    anglel = np.radians(anglel_deg)
    angler = np.radians(angler_deg)
    randoma = np.radians(randomangle)

    # Generation 0: single vertical trunk
    generations_data = [
        [(0.0, 0.0, initial_length, np.pi / 2)]
    ]

    for g in range(1, generations + 1):
        new_generation = []

        for x, y, length, theta in generations_data[g - 1]:
            # End point of current branch
            x_end = x + length * np.cos(theta)
            y_end = y + length * np.sin(theta)

            # Left branch
            new_length = ((1-randombranch)*length * scale 
                          + randombranch * length * scale * (1-2*np.random.rand()) )
                   
            angleleft = anglel + randoma*(1-2*np.random.rand())
            new_generation.append(
                (x_end, y_end, new_length, theta + angleleft)
            )
           
            # Right branch            
            new_length = ((1-randombranch)*length * scale 
                          + randombranch * length * scale * (1-2*np.random.rand()) )        

            angleright = angler + randoma*(1-2*np.random.rand())
            new_generation.append(
                (x_end, y_end, new_length, theta - angleright)
            )

        generations_data.append(new_generation)

    return generations_data

# ------------------------------
#
# Generate Tree
#
# ------------------------------
tree = generate_sierpinski_tree(
    generations=generations,
    initial_length=1.0,
    anglel_deg=angleleft,
    angler_deg=angleright,
    scale=branchscale
)

# ---------------------------
# Define Plot
# ---------------------------
def gradient_colors(n, cmap_name="viridis"):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(i) for i in np.linspace(0, 1, n)]

colormap = gradient_colors(generations+1, "plasma")


fig, ax = plt.subplots(1, 1, figsize=(5, 5),
                        layout="constrained")
ax.set_aspect('equal')
ax.axis('off')
ax.set(xlim=(-2,2),ylim=(0,4))

# info box
infobox = name+'\n'
infobox +='generations: ' + "{:.0f}".format(generations)+'\n'
infobox +='angle left: ' + "{:.0f}".format(angleleft)+'\n'
infobox +='angle right: ' + "{:.0f}".format(angleright)+'\n'
infobox +='branch scale: ' + "{:.2f}".format(branchscale)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)

x0 = 0  
y0 = 2
xscale = 2
yscale = 2

# -----------------------------------------------------------------------------
#
# Main Routine
#
# -----------------------------------------------------------------------------
def animate(k):
    
  kindex = int(k/sweep)  
  print(kindex,' von ',generations,' branches ',len(tree[kindex])) 
  
  ax.clear()
  ax.set_aspect('equal')
  ax.axis('off')
  move = k/(generations*sweep)
  x0 = 0 - move*0.5 
  y0 = 2
  xscale = 1+1*(1-move)
  yscale = 1+1*(1-move)
  xmin = x0-xscale
  xmax = x0+xscale
  ymin = y0-yscale
  ymax = y0+yscale 
  ax.set(xlim=(xmin,xmax),ylim=(ymin,ymax))

  ax.text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)

  if kindex>0:
    for i in range(kindex):
      generation = tree[i]    
      for x, y, length, theta in generation:
        ll = length
        x_end = x + ll * np.cos(theta)
        y_end = y + ll * np.sin(theta)
        ax.plot([x, x_end], [y, y_end], color=colormap[i], linewidth=2-1*i/generations)
      
  generation = tree[kindex]    
  for x, y, length, theta in generation:
        ll = length*(k-kindex*sweep)/(sweep)
        x_end = x + ll * np.cos(theta)
        y_end = y + ll * np.sin(theta)
        ax.plot([x, x_end], [y, y_end], color=colormap[kindex], linewidth=2-1*kindex/generations,alpha=(k-kindex*sweep)/(sweep))

  fig.suptitle('Sierpinski Tree, Generation: '+"{:.0f}".format(kindex))


anim = animation.FuncAnimation(fig,animate,interval=1,frames=generations*sweep)
anim.save(animationname,fps=25,dpi=300)       
    
    
