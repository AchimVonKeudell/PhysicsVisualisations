#
# Cross section
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

model = "Coulomb"
#model = "Langevin"

el = 1.602176634e-19 # el charge
amu = 1.6605390666e-27 # proton mass
kB = 1.380649e-23 # Boltzmann const.
eps0 = 8.854e-12 # epsilon 0
me = 9e-31 # mass electron

# Initialize Simulation
if model == "Langevin":
    m = 4
    animationname = 'cross_langevin.mp4'
if model == "Coulomb":
    m = 1
    animationname = 'cross_coulomb.mp4'    

t = 0 
dt = 1e-11 
dtplot = 5e-10
tend = 30*74e-9
steps = int(tend/dtplot)      

tspan =[0]
Espan = [0]

# setup plot
xmax = 10
ymax = 10

fig, ax = plt.subplots(1,2,figsize=(8,4.5))
#plt.subplots_adjust(left=0.1, right=0.9, top=0.65, bottom=0.1)

linea1 = []
linea2 = []
linea3 = []
for i in range(100):
  linea1.append(ax[0].plot(0,0,color='b',lw=1, label ='')[0])
  linea2.append(ax[0].plot(0,0,color='g',lw=1, label ='')[0])
  linea3.append(ax[0].plot(0,0,color='r',lw=1, label ='')[0])
ax[0].set(xlabel="x",ylabel="y")
ax[0].set(xlim=(-xmax,xmax),ylim=(-ymax,ymax))
lineat = ax[0].plot(0,0,color='r',lw=0, label ='',marker='o')[0]
linea1[0].set_label('E0')
linea2[0].set_label('E0 * 0.5')
linea3[0].set_label('E0 * 0.25')
ax[0].legend(fontsize = 6,loc=1)

# setup info box left
infobox = ''
infobox += 'V ~ 1/r^n, n: ' + "{:.0f}".format(m) +'\n'
infobox += 'V nature: attractive'
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[0].text(0.05,0.1,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax[0].transAxes)


liner1 = []
liner2 = []
liner3 = []
for i in range(100):
  liner1.append(ax[1].plot(0,0,color='b',lw=1, label ='')[0])
  liner2.append(ax[1].plot(0,0,color='g',lw=1, label ='')[0])
  liner3.append(ax[1].plot(0,0,color='r',lw=1, label ='')[0])
ax[1].set(xlabel="x",ylabel="y")
ax[1].set(xlim=(-xmax,xmax),ylim=(-ymax,ymax))
linert = ax[1].plot(0,0,color='r',lw=0, label ='',marker='o')[0]
liner1[0].set_label('E0')
liner2[0].set_label('E0 * 0.5')
liner3[0].set_label('E0 * 0.25')
ax[1].legend(fontsize = 6,loc=1)


# setup info box left
infobox = ''
infobox += 'V ~ 1/r^n, n: ' + "{:.0f}".format(m)+'\n' 
infobox += 'V nature: repulsive' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[1].text(0.05,0.1,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax[1].transAxes)

fig.tight_layout(pad=1,rect=(0,0,1,0.95))

dt = 0.0001
dtreplot = 0.001
Ekin = 0.00001
 
def force(r):
      rb = np.linalg.norm(r)
      return -1/rb**(m+1)*1/rb*r     

def animate(k):
    
    global IEDF,Escan,Particles,IEDFNorm,t,a,x,v,ax
    
    # -----------------------------------------------------------------------
    # attractive potential
    # -----------------------------------------------------------------------
        
    # generate trajectory 1
    r1 = np.array([-10,8-k/49*8,0])
    v0 = np.sqrt(2*Ekin*el/amu)
    v1 = np.array([v0,0,0])
    a1 = np.array([0,0,0])
    t = 0
    time = 0
    xt = []
    yt = []
    while r1[0]>=-xmax and r1[0]<=xmax and r1[1]>=-ymax and r1[1]<=ymax and t<30:
        t += dt
        time += dt      
        dr1 = dt*v1+0.5*dt**2*a1
        r1 = r1+dr1
        v1 = v1 + 0.5*dt*a1
        a1 = force(r1)*1000
        if np.linalg.norm(a1)>1000000: break
        v1 = v1 + 0.5*dt*a1      
        if np.linalg.norm(r1)<0.1: break
        if time>dtreplot:
            time = 0            
            xt.append(r1[0])
            yt.append(r1[1])
    linea1[k].set_xdata(xt)
    linea1[k].set_ydata(yt) 
    
    # generate trajectory 2
    r1 = np.array([-10,8-k/49*8,0])
    v0 = np.sqrt(2*Ekin*el/amu)*0.5
    v1 = np.array([v0,0,0])
    a1 = np.array([0,0,0])
    t = 0
    time = 0
    xt = []
    yt = []
    while r1[0]>=-xmax and r1[0]<=xmax and r1[1]>=-ymax and r1[1]<=ymax and t<30:
        t += dt
        time += dt      
        dr1 = dt*v1+0.5*dt**2*a1
        r1 = r1+dr1
        v1 = v1 + 0.5*dt*a1
        a1 = force(r1)*1000
        if np.linalg.norm(a1)>1000000: break
        v1 = v1 + 0.5*dt*a1      
        if np.linalg.norm(r1)<0.1: break
        if time>dtreplot:
            time = 0            
            xt.append(r1[0])
            yt.append(r1[1])
    linea2[k].set_xdata(xt)
    linea2[k].set_ydata(yt) 
    
    # generate trajectory 3
    r1 = np.array([-10,8-k/49*8,0])
    v0 = np.sqrt(2*Ekin*el/amu)*0.25
    v1 = np.array([v0,0,0])
    a1 = np.array([0,0,0])
    t = 0
    time = 0
    xt = []
    yt = []
    while r1[0]>=-xmax and r1[0]<=xmax and r1[1]>=-ymax and r1[1]<=ymax and t<30:
        t += dt
        time += dt      
        dr1 = dt*v1+0.5*dt**2*a1
        r1 = r1+dr1
        v1 = v1 + 0.5*dt*a1
        a1 = force(r1)*1000
        if np.linalg.norm(a1)>1000000: break
        v1 = v1 + 0.5*dt*a1      
        if np.linalg.norm(r1)<0.1: break
        if time>dtreplot:
            time = 0            
            xt.append(r1[0])
            yt.append(r1[1])
    linea3[k].set_xdata(xt)
    linea3[k].set_ydata(yt) 
    
    
    # -----------------------------------------------------------------------
    # repulsive potential
    # -----------------------------------------------------------------------
    
    # generate trajectory 1
    r1 = np.array([-10,8-k/49*8,0])
    v0 = np.sqrt(2*Ekin*el/amu)
    v1 = np.array([v0,0,0])
    a1 = np.array([0,0,0])
    t = 0
    time = 0
    xt2 = []
    yt2 = []
    while r1[0]>=-xmax and r1[0]<=xmax and r1[1]>=-ymax and r1[1]<=ymax and t<30:
        t += dt
        time += dt
        dr1 = dt*v1+0.5*dt**2*a1
        r1 = r1+dr1
        v1 = v1 + 0.5*dt*a1
        a1 = -force(r1)*1000
        v1 = v1 + 0.5*dt*a1
        if time>dtreplot:
            time = 0
            xt2.append(r1[0])
            yt2.append(r1[1])
    liner1[k].set_xdata(xt2)
    liner1[k].set_ydata(yt2) 
    
    # generate trajectory 2
    r1 = np.array([-10,8-k/49*8,0])
    v0 = np.sqrt(2*Ekin*el/amu)*0.5
    v1 = np.array([v0,0,0])
    a1 = np.array([0,0,0])
    t = 0
    time = 0
    xt2 = []
    yt2 = []
    while r1[0]>=-xmax and r1[0]<=xmax and r1[1]>=-ymax and r1[1]<=ymax and t<30:
        t += dt
        time += dt
        dr1 = dt*v1+0.5*dt**2*a1
        r1 = r1+dr1
        v1 = v1 + 0.5*dt*a1
        a1 = -force(r1)*1000
        v1 = v1 + 0.5*dt*a1
        if time>dtreplot:
            time = 0
            xt2.append(r1[0])
            yt2.append(r1[1])
    liner2[k].set_xdata(xt2)
    liner2[k].set_ydata(yt2) 
    
    # generate trajectory 3
    r1 = np.array([-10,8-k/49*8,0])
    v0 = np.sqrt(2*Ekin*el/amu)*0.25
    v1 = np.array([v0,0,0])
    a1 = np.array([0,0,0])
    t = 0
    time = 0
    xt2 = []
    yt2 = []
    while r1[0]>=-xmax and r1[0]<=xmax and r1[1]>=-ymax and r1[1]<=ymax and t<30:
        t += dt
        time += dt
        dr1 = dt*v1+0.5*dt**2*a1
        r1 = r1+dr1
        v1 = v1 + 0.5*dt*a1
        a1 = -force(r1)*1000
        v1 = v1 + 0.5*dt*a1
        if time>dtreplot:
            time = 0
            xt2.append(r1[0])
            yt2.append(r1[1])
    liner3[k].set_xdata(xt2)
    liner3[k].set_ydata(yt2) 
    
    print('Collision parameter: '+"{:.2f}".format(5-k/100*4.5))
    fig.suptitle('Collision parameter: '+"{:.2f}".format(5-k/100*4.5))
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=50)
anim.save(animationname,fps=5,dpi=180)