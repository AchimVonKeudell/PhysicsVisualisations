#
# Mandelbrot Fractal
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2026
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation  as animation



mode = 'Tip'
#mode = 'Bottom Branch'
#mode = 'Seahorse'

if mode == 'Tip':
  name = 'Tip Point'
  x0 = -2
  y0 = 0
  modus = 'animation'
  #modus = 'endpoint'
  xscale0 =3
  yscale0 = 3
  iterations = 300
  ordersofmagnitude = 10
  maxiterations = 400
  animationname = 'mandelbrot_tip.gif'

if mode == 'Bottom Branch':
  name = 'Bottom Point'
  x0 = 0
  y0 = 1
  modus = 'animation'
  #modus = 'endpoint'
  xscale0 =3
  yscale0 = 3
  iterations = 300
  ordersofmagnitude = 10
  maxiterations = 400
  animationname = 'mandelbrot_bottom.gif'

if mode == 'Seahorse':
  name = 'Seahorse Point'
  x0 = -0.7494
  y0 = 0.1004
  modus = 'animation'
  #modus = 'endpoint'
  xscale0 =3
  yscale0 = 3
  iterations = 300
  ordersofmagnitude = 10
  maxiterations = 400
  animationname = 'mandelbrot_seahorse.gif'


def mandelbrot(width=800, height=800, max_iter=100,xscale=1,yscale=1):
    # Define region in complex plane
    re = np.linspace(x0-xscale, x0+xscale, width)
    im = np.linspace(y0-yscale, y0+yscale, height)
    c = re[np.newaxis, :] + 1j * im[:, np.newaxis]

    z = np.zeros_like(c)
    img = np.zeros(c.shape, dtype=int)

    for i in range(max_iter):
        mask = np.abs(z) <= 2
        z[mask] = z[mask]**2 + c[mask]
        img[mask] = i

    return img

# Generate Mandelbrot set

fig, ax = plt.subplots(1, 1, figsize=(5, 5),
                        layout="constrained")


ax.set_aspect('equal')
ax.axis('off')

# info box
infobox = name + ' \n'
infobox +='x$_0$: ' + "{:.0f}".format(x0) + '\n '
infobox +='y$_0$: ' + "{:.0f}".format(y0)

props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8) 
ax.text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)
  
    

def animate(k):
      
  print(k, 'of ',iterations)  
  xscale = xscale0*10**(-ordersofmagnitude*k/iterations)
  yscale = yscale0*10**(-ordersofmagnitude*k/iterations)
  
  iterat = 100 + k/iterations*(maxiterations-100)
  mandelbrot_set = mandelbrot(max_iter=maxiterations,xscale=xscale,yscale=yscale)

  ax.imshow(np.flip(mandelbrot_set,0), cmap="plasma", extent=[x0-xscale, x0+xscale, y0-yscale, y0+yscale])
  fig.suptitle('Mandelbrot set, Magnification 10$^x$  x: '+"{:.2f}".format(ordersofmagnitude*k/iterations))

if modus == 'endpoint':
  xscale = xscale0*10**(-ordersofmagnitude)
  yscale = yscale0*10**(-ordersofmagnitude)
  
  mandelbrot_set = mandelbrot(max_iter=maxiterations,xscale=xscale,yscale=yscale)
  ax.imshow(np.flip(mandelbrot_set,0), cmap="plasma", extent=[x0-xscale, x0+xscale, y0-yscale, y0+yscale])
    
else:    
    
  anim = animation.FuncAnimation(fig,animate,interval=1,frames=iterations)
  anim.save(animationname,fps=25,dpi=300)   
