#
# Cavitation
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
import matplotlib.animation as animation

model = 'liquid'

if model == 'liquid':
  tend = 2e-4
  steps = 400
  R0 = 25e-6 # initial radius HV tip 
  pressure0 =  1e9 # pressure in the plasma after ignition
  UkV0 = 20000 # voltage HV pin in V, used to cacluate electric field pressure
  tpulse = 10e-9 # pulse length in s
  vshock = 70 # acoustic shockwave in m/s
  animationname = 'bubble.mp4'


# parameter wasser 
sigmawater = 0.072  # surface tension in Nm^-2 
etawater = 0.001  # visocsity in Pa s 
pressureeq = 3e5 # pressure infty in Pa 
rhoeq = 1000 # density liquid in kg m^-3 
gamma = 1.95 # adiabatic exponent gas 
epsilon0 = 8.854e-12 # dielectric constant 
Bwater = 3000*1e5 # equation of state water 
nwater = 7 # equation of state water 
cwater0 = 1435 # sound velocity in ms^-1 
kB = 1.38e-23 # Boltzmann constant

# start conditions 
U0 = 0 # initial velocity 

#pplasma = 3/(1.6e-19)*3e28;
T0 = pressure0/(3e28*kB) # assuming water density at the beginning
print('Initial Temperature assuming liquid water: ',int(T0),' (K)')

# Efield and Efield induced pressure
def Efield(UkV):
    return UkV/R0

def pEfield(UkV,t):
    if t <tpulse:
        return 0.5*epsilon0*Efield(UkV)**2
    else:
        return 0
     
# pressure at the interface
def pgas(R,Rdot):
    return pressure0*(R0/R)**(3*gamma) - 2*sigmawater/R - 4*etawater/R*Rdot

# sound velcoiyt water 
Bwater = 3000*1e5
nwater = 7
def cwater(R,Rdot): 
  return cwater0*((pgas(R, Rdot) + Bwater)/(Bwater + pressureeq))**((nwater - 1)/(2*nwater))

# enthalpy 
def enthalpy(R,Rdot,t): 
  return (1/rhoeq*nwater/(nwater - 1)*(Bwater + 
     pressureeq)*(((pEfield(UkV0,t) + pgas(R,Rdot) + 
          Bwater)/(Bwater + pressureeq))**((nwater - 1)/nwater) - 1))

# solving the RP equation 
def func(y,t):
    R = y[0]
    Rp = y[1]
    Rpp = 1/R*1/(1-Rp/cwater(R,Rp))*(-1.5*(1-Rp/(3*cwater(R,Rp)))*Rp**2+(1+Rp/cwater(R,Rp))*enthalpy(R,Rp,t))
    return [Rp,Rpp]

# Mathematica version NDsolve
#lsg = NDSolve[{(1 - R'[t]/cwater[R[t], R'[t]]) R[t] R''[t] + 
#     1.5 (1 - R'[t]/(3 cwater[R[t], R'[t]])) (R'[t])^2 == (1 + 
#       R'[t]/cwater[R[t], R'[t]]) enthalpy[R[t], R'[t], t], 
#   R'[0] == U0, R[0] == R0}, R, {t, 10^-10, 10^-3}];
#T[t_] = T0 (R0/R[t])^(3 (gamma - 1));

# Calculate solution
print('Solve Rayleigh-Plesset equation')
y0 = [R0,U0]
t = np.linspace(1e-10,tend,steps)
y = odeint(func,y0,t)

# -----------------------------------------------------------------------
# Define Plot
# -----------------------------------------------------------------------
# dynamic
fig, ax = plt.subplots(1,2,figsize=(8,4.5))
tp = []
r = []
temp = []
line = ax[0].plot(tp,r,label='Radius')[0]
axp = ax[0].twinx()
axp.set_yscale('log')
linetemp = axp.plot(tp,temp,color='orange',label='Temperature')[0]
axp.set(ylabel='T (K)',ylim=(1,5000))
ax[0].set(xlabel='t (s)',ylabel='R (micrometer)',xlim=(0,tend/1e-6),ylim=(0,500))
ax[0].legend(fontsize=6)
axp.legend(fontsize=6,loc=4)

# illustration
ax[1].set_aspect(1)
ax[1].set(xlabel='x (micrometer)',ylabel='y (micrometer)',xlim=(-500,500),ylim=(-500,500))  
c1 = ax[1].add_patch(Rectangle((-500,-500),1000,1000,color='lightblue'))
c2 = ax[1].add_patch(Circle((0,0),25e-6/1e-6,color=[1,1,0]))
textt = ax[1].text(0.5,0.5,'T', fontsize=8,color='w',horizontalalignment='center',verticalalignment='center',transform=ax[1].transAxes)

# shockwave
circlex = []
circley = []
angle = np.linspace(0,2*np.pi,360)
circlex = R0/1e-6*np.cos(angle)
circley = R0/1e-6*np.sin(angle)
linecircle = ax[1].plot(circlex,circley,'g--',label='shockwave - 0.01 v_sound')[0]
ax[1].legend(fontsize=6,loc=1)

# info box
infobox = ''
infobox += 'p0: ' + "{:.0e}".format(pressure0) + ' (Pa)\n' 
infobox += 'R0: ' + "{:.0f}".format(R0/1e-6) + ' (mu)\n' 
infobox += 'T0: ' + "{:.0f}".format(T0) + ' (K)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[0].text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax[0].transAxes)

fig.tight_layout(pad=3)

# ----------------------------------
# Main loop
# ----------------------------------
dt = tend/steps
def animate(k):
   global c1,c2 
   tp.append(t[k]/1e-6)
   r.append(y[k,0]/1e-6)
   line.set_xdata(tp)
   line.set_ydata(r)
   
   T = T0*(R0/y[k,0])**(3*(gamma - 1))
   temp.append(T)
   linetemp.set_xdata(tp)
   linetemp.set_ydata(temp)
 
   #circlex = y[k,0]/1e-6*np.cos(angle)
   #circley = y[k,0]/1e-6*np.sin(angle)
   circlex = vshock*t[k]/1e-6*np.cos(angle)
   circley = vshock*t[k]/1e-6*np.sin(angle)
   linecircle.set_xdata(circlex)
   linecircle.set_ydata(circley)
   
   colort = (T/T0)**0.05 # choose a nonlinear color scale
   #print(colort)
   if colort>1: colort = 1   
   c1.remove()
   c2.remove()   
   c1 = ax[1].add_patch(Rectangle((-500,-500),1000,1000,color='lightblue'))
   c2 = ax[1].add_patch(Circle((0,0),y[k,0]/1e-6,color=([colort,colort,0])))
   textt.set_text('T: ' + "{:.1f}".format(T)+ ' K')

   fig.suptitle('time: ' + "{:.0f}".format(t[k]/1e-6)+' microseconds')
   print('Animate frame ',k,' of ',steps)
    
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)