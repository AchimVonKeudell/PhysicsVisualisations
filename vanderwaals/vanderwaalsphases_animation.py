#
# van der Waals gas phases
# 
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#
# initial source code from https://scipython.com/blog/the-maxwell-construction/
#
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.optimize import newton
from scipy.signal import argrelextrema
import matplotlib.animation as animation
from matplotlib.patches import Rectangle



gas = 'CO$_2$'
# Critical pressure, volume and temperature
# These values are for the van der Waals equation of state for CO2:
# (p - a/V^2)(V-b) = RT. Units: p is in Pa, Vc in m3/mol and T in K.
pc = 7.404e6
Vc = 1.28e-4
Tc = 304

Tr = 270/Tc


def vdw(Tr, Vr):
    """Van der Waals equation of state.

    Return the reduced pressure from the reduced temperature and volume.

    """

    pr = 8*Tr/(3*Vr-1) - 3/Vr**2
    return pr

def get_Vlims_R(pr0):
    """Solve the inverted van der Waals equation for reduced volume.

        Return the lowest and highest reduced volumes such that the reduced
        pressure is pr0. It only makes sense to call this function for
        T<Tc, ie below the critical temperature where there are three roots.

    """

    eos = np.poly1d( (3*pr0, -(pr0+8*Tr), 9, -3) )
    roots = eos.r
    roots.sort()
    Vrmin, _, Vrmax = roots
    return Vrmin, Vrmax


def vdw_maxwell(Tr, Vr):
    """Van der Waals equation of state with Maxwell construction.

    Return the reduced pressure from the reduced temperature and volume,
    applying the Maxwell construction correction to the unphysical region
    if necessary.

    """

    pr = vdw(Tr, Vr)
    if Tr >= 1:
        # No unphysical region above the critical temperature.
        return pr

    if min(pr) < 0:
         raise ValueError('Negative pressure results from van der Waals'
                         ' equation of state with Tr = {} K.'.format(Tr))

    # Initial guess for the position of the Maxwell construction line:
    # the volume corresponding to the mean pressure between the minimum and
    # maximum in reduced pressure, pr.
    iprmin = argrelextrema(pr, np.less)
    iprmax = argrelextrema(pr, np.greater)
    Vr0 = np.mean([Vr[iprmin], Vr[iprmax]])

    def get_Vlims(pr0):
        """Solve the inverted van der Waals equation for reduced volume.#

        Return the lowest and highest reduced volumes such that the reduced
        pressure is pr0. It only makes sense to call this function for
        T<Tc, ie below the critical temperature where there are three roots.

        """

        eos = np.poly1d( (3*pr0, -(pr0+8*Tr), 9, -3) )
        roots = eos.r
        roots.sort()
        Vrmin, _, Vrmax = roots
        return Vrmin, Vrmax

    def get_area_difference(Vr0):
        """Return the difference in areas of the van der Waals loops.

        Return the difference between the areas of the loops from Vr0 to Vrmax
        and from Vrmin to Vr0 where the reduced pressure from the van der Waals
        equation is the same at Vrmin, Vr0 and Vrmax. This difference is zero
        when the straight line joining Vrmin and Vrmax at pr0 is the Maxwell
        construction.

        """

        pr0 = vdw(Tr, Vr0)
        Vrmin, Vrmax = get_Vlims(pr0)
        return quad(lambda vr: vdw(Tr, vr) - pr0, Vrmin, Vrmax)[0]

    # Root finding by Newton's method determines Vr0 corresponding to
    # equal loop areas for the Maxwell construction.
    Vr0 = newton(get_area_difference, Vr0)
    pr0 = vdw(Tr, Vr0)
    Vrmin, Vrmax = get_Vlims(pr0)

    # Set the pressure in the Maxwell construction region to constant pr0.
    pr[(Vr >= Vrmin) & (Vr <= Vrmax)] = pr0
    return pr, Vrmin, Vrmax

vsteps = 500
Vr = np.linspace(0.5, 4, vsteps)
pr = np.ones(vsteps)
prz = np.zeros(vsteps)

fig, ax = plt.subplots(1,2,figsize=(9,4.5))
fig.suptitle('van der Waals gas')
ax[0].set_xlim(0, 4)
ax[0].set_xlabel('V/V$_r$')
ax[0].set_ylim(0, 1.6)
ax[0].set_ylabel('p/p$_r$')

areax1 = [1]
areay1 = [1]
areax2 = [1]
areay2 = [1]

textT = ax[0].text(0.7,0.74,'T$_r$ = {:.2f}'.format(Tr),transform=ax[0].transAxes)
textT1 = ax[0].text(0.7,0.68,'T = {:.0f}'.format(Tr*Tc)+' K',transform=ax[0].transAxes)

textfl = ax[1].text(1,1,'')
textg = ax[1].text(1,1,'')

ax[1].axis('off')
ax[1].set(ylim=(0,0.6),xlim=(-5,5))
patchl = ax[1].add_patch(Rectangle((-1,0),2,0.5,ec='b',fc='blue'))
patchg = ax[1].add_patch(Rectangle((-1,0),2,0.5,ec='b',fc='lightblue'))
piston = ax[1].add_patch(Rectangle((-1,0),2,0.01,ec='black',fc='black'))
piston1 = ax[1].add_patch(Rectangle((-0.05,0),0.1,0.1,ec='black',fc='black'))
   
ax[1].add_patch(Rectangle((-1.2,0),0.2,0.6,ec='lightgrey',fc='grey'))
ax[1].add_patch(Rectangle((1,0),0.2,0.6,ec='lightgrey',fc='grey'))

textfl = ax[1].text(1,1,'')
textg = ax[1].text(1,1,'')


# info box
infobox = ''
infobox += 'gas: ' + gas + '\n' 
infobox += 'p$_c$: ' + "{:.2f}".format(pc/1e5) + ' (bar)\n' 
infobox += 'T$_c$: ' + "{:.0f}".format(Tc) + ' (K)\n' 
infobox += 'V$_c$: ' + "{:.0f}".format(Vc*1e6) + ' (cm$^3$/mol)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[0].text(0.7,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax[0].transAxes)

for i in range(1000):

     Tup = 350
     Tdown = 257
     Tr = 1 - (i+1)/1000*(Tc-Tdown)/Tc
    
     Vmin = vdw_maxwell(Tr, Vr)[1]
     pmin = vdw(Tr,Vmin)
     areax1.append(Vmin)
     areay1.append(pmin)
     
     Vmax = vdw_maxwell(Tr, Vr)[2]
     pmax = vdw(Tr,Vmax)
     areax2.append(Vmax)
     areay2.append(pmax)
     

line4 = ax[0].fill_between(areax1,areay1,color='orange',alpha=0.3,edgecolor=None)  
line5 = ax[0].fill_between(areax2,areay2,color='orange',alpha=0.3,edgecolor=None)

line1 = ax[0].plot(Vr, vdw(270/Tc,Vr), lw=1, color='red',linestyle='dashed')[0]
line2 = ax[0].plot(Vr, vdw_maxwell(270/Tc, Vr)[0], lw=1, color='red')[0]
line3 = ax[0].fill_between(Vr,vdw(270/Tc,Vr),vdw_maxwell(270/Tc, Vr)[0],color='red',lw=0,alpha=0.3)
 
line6 = ax[0].plot([],[],marker='o',color='black',markersize=5,lw=0)[0]


pr = vdw_maxwell(270/Tc, Vr)[0]
Vmin = vdw_maxwell(270/Tc, Vr)[1]
Vmax = vdw_maxwell(270/Tc, Vr)[2]

def animate(k):
   global line3,line4,line5
   global line1,line2
   global textflg
   global textfl
   global textg
   global patchl,patchg,piston,piston1
 
   Vred = Vr[500-k-1]*0.5/4
   
   patchl.remove()
   patchg.remove()
   piston.remove()
   piston1.remove()
   if Vr[500-k-1]>Vmax:
     patchl = ax[1].add_patch(Rectangle((-1,0),2,0,ec=None,fc='blue'))
     patchg = ax[1].add_patch(Rectangle((-1,0),2,Vred,ec=None,fc='lightblue'))
     
     textg.remove()
     textg = ax[1].text(0,Vred/2,'gas',horizontalalignment='center',verticalalignment='center',fontsize=8)
   
   if Vmin<Vr[500-k-1]<Vmax:
     fg = (Vr[500-k-1]-Vmin)/(Vmax-Vmin)
     fl = 1-fg
     
     patchl = ax[1].add_patch(Rectangle((-1,0),2,Vred*fl,ec=None,fc='blue'))
     patchg = ax[1].add_patch(Rectangle((-1,Vred*fl),2,Vred*fg,ec=None,fc='lightblue'))
     
     textg.remove()
     textg = ax[1].text(0,Vred*fl+(Vred*fg)/2,'gas',horizontalalignment='center',verticalalignment='center',fontsize=8)   
     textfl.remove()
     textfl = ax[1].text(0,Vred*fl/2,'liquid',horizontalalignment='center',verticalalignment='center',fontsize=8,color='white')
 
   if Vr[500-k-1]<Vmin:
     patchl = ax[1].add_patch(Rectangle((-1,0),2,Vred,ec=None,fc='blue'))
     patchg = ax[1].add_patch(Rectangle((-1,0),2,0,ec=None,fc='lightblue'))
     
     textg.remove()
     textfl.remove()
     textfl = ax[1].text(0,Vred/2,'liquid',horizontalalignment='center',verticalalignment='center',fontsize=8,color='white')
 
   piston = ax[1].add_patch(Rectangle((-1,Vred),2,0.01,ec='black',fc='black'))
   piston1 = ax[1].add_patch(Rectangle((-0.05,Vred),0.1,0.1,ec='black',fc='black'))

   
   line6.set_xdata([Vr[500-k-1]])
   line6.set_ydata([pr[500-k-1]])
  
anim = animation.FuncAnimation(fig,animate,interval=1,frames=vsteps)
anim.save('vanderWaalsphases.mp4',fps=25,dpi=300)