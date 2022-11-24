#
# Fermi acceleartion of an electron
# bouncing off a matrix sheath
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

model = "NoCollisions"
#model = "Collsions"

el = 1.602176634e-19 # el charge
amu = 1.6605390666e-27 # proton mass
kB = 1.380649e-23 # Boltzmann const.
eps0 = 8.854e-12 # epsilon 0
me = 9e-31 # mass electron
Te = 3 # electron temperature

# Initialize Simulation
if model == "NoCollisions":
    nu = 0 # collision frequency 
    Emax = 30
if model == "Collisions":
    nu = 1e7 # collision frequency        
    Emax = 30
V0 = 700 # Peak to peak
fRF = 13.6 # in MHz
omegaRF = 2*np.pi*fRF*1e6 # angular frequency
zwidth = 0.04 # width of simulation span in m
s0 = 0.004 # average sheath width
n0 = V0*eps0/el*2*1/(2*s0)**2 # charge density in space charge sheath
v = [np.sqrt(Te*el/(me)),0,0] # initial velocity thermal

x = 0
y = 0
t = 0 
dt = 1e-11 
dtplot = 5e-10
tend = 30*74e-9
steps = int(tend/dtplot)      

tspan =[0]
Espan = [0]

# setup plot
sleftx = [-zwidth/2/1e-3,-zwidth/2/1e-3]
slefty = [-1,1]
srightx = [zwidth/2/1e-3,zwidth/2/1e-3]
srighty = [-1,1]

fig, ax = plt.subplots(1,2,figsize=(8,4.5))
plt.subplots_adjust(left=0.1, right=0.8, top=0.85, bottom=0.1)
line1 = ax[0].plot(x,y,color='b',lw=1, label ='',marker='o')[0]
lineleft = ax[0].plot(sleftx,slefty,color='r',lw=1)[0] 
lineright = ax[0].plot(srightx,srighty,color='r',lw=1)[0] 
ax[0].set(xlabel="x (mm)",ylabel="y (mm)")
ax[0].set(xlim=(-zwidth/2/1e-3,+zwidth/2/1e-3),ylim=(-1,1))

lineE = ax[1].plot(tspan,Espan,color='g',lw=1, label ='')[0]
ax[1].set(xlabel="time (ns)",ylabel="energy (eV)")
ax[1].set(xlim=(0,tend/1e-9),ylim=(0,Emax))
fig.tight_layout(pad=2)

# setup info box
infobox = ''
infobox += 'V0: ' + "{:.0f}".format(V0) + ' (V)\n' 
infobox += 'fRF: ' + "{:.2f}".format(fRF) + ' (Hz)\n' 
infobox += 'n0: ' + "{:.2e}".format(n0) + ' (m^-3)\n' 
infobox += 'nu: ' + "{:.0e}".format(nu) + ' (1/s)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[0].text(0.05,0.95,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax[0].transAxes)

# -----------------------------------------------------------------------------
# Collision calculation electrons
#  from eduPIC code example by Z. Donko et al.
# -----------------------------------------------------------------------------
def collision(v):
        vb = np.sqrt(v[0]**2+v[1]**2+v[2]**2)
        # calulate Euler angles
        if v[0] == 0:
            theta = 0.5*np.pi
        else:
            theta = np.arctan(np.sqrt(v[1]**2+v[2]**2)/v[0])
        if v[1] == 0:
            if v[2] > 0:
                phi = 0.5*np.pi
            else:
                phi = -0.5*np.pi
        else:        
            phi = np.arctan(v[2]/v[1])

        st = np.sin(theta)
        ct = np.cos(theta)
        sp = np.sin(phi)
        cp = np.cos(phi)

        # select scattering angles 
        xi = np.arccos(1-2*np.random.rand())
        eta = 2*np.pi*np.random.rand()
        
        sc = np.sin(xi)
        cc = np.cos(xi)
        se = np.sin(eta)
        ce = np.cos(eta)
        
        # rotate velocity vector
        vx = vb * (ct * cc - st * sc * ce);
        vy = vb * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
        vz = vb * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
   
        return np.array([vx,vy,vz]) 

# -----------------------------------------------------------------------------
# Calculate E field in matrix sheath
# -----------------------------------------------------------------------------
def Efeld(x,t):
    # position sheath edges
    s1 = -zwidth/2+s0*(1-np.cos(omegaRF*t))
    s2 = +zwidth/2-s0*(1+np.cos(omegaRF*t))
    
    # E field Matrix sheath
    if x<s1:
        return -0.5*el/eps0*n0*abs(x-s1)
    elif x>s2:
        return 0.5*el/eps0*n0*abs(x-s2)
    else:
        return 0

def animate(k):
    
    global IEDF,Escan,Particles,IEDFNorm,t,a,x,v
    
    time = 0
    while time<dtplot:
        t += dt
        time += dt
        a = -el*Efeld(x,t)/me   
        v[0] += a*dt # propagate velocity
         # check after microsteps*dt whether a collsion should take place
        if np.random.rand() < 1-np.exp(-nu*dt):
           v = collision(v)
        x += v[0]*dt # new position        
    
    # position sheath edges
    s1 = -zwidth/2+s0*(1-np.cos(omegaRF*t))
    s2 = +zwidth/2-s0*(1+np.cos(omegaRF*t))
    
    sleftx = [s1/1e-3,s1/1e-3]      
    srightx = [s2/1e-3,s2/1e-3]      
    line1.set_xdata(x/1e-3)
    line1.set_ydata(y) 
    lineleft.set_xdata(sleftx)
    lineleft.set_ydata(slefty)
    lineright.set_xdata(srightx)
    lineright.set_ydata(srighty)
    
    Energy = 0.5*me*(v[0]**2+v[1]**2+v[2]**2)/el
    tspan.append(t/1e-9)
    Espan.append(Energy)
    lineE.set_xdata(tspan)
    lineE.set_ydata(Espan) 
    
    print('time: '+"{:.1f}".format(t/1e-9)+' ns of '+"{:.0f}".format(tend/1e-9))
    fig.suptitle('time: ' + "{:.0f}".format(t/1e-9)+' ns')
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save('fermi.gif',fps=25,dpi=180)