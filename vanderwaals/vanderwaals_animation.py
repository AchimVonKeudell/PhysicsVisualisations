#
# van der Waals gas
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


gas = 'CO$_2$'
# Critical pressure, volume and temperature
# These values are for the van der Waals equation of state for CO2:
# (p - a/V^2)(V-b) = RT. Units: p is in Pa, Vc in m3/mol and T in K.
pc = 7.404e6
Vc = 1.28e-4
Tc = 304

Tr = 0.9


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

Vr = np.linspace(0.5, 4, 500)
pr = np.ones(500)
prz = np.zeros(500)

fig, ax = plt.subplots(figsize=(4.5,4.5))
fig.suptitle('van der Waals gas')
ax.set_xlim(0, 4)
ax.set_xlabel('V/V$_r$')
ax.set_ylim(0, 1.6)
ax.set_ylabel('p/p$_r$')
#ax.legend(fontsize=6)

line1 = ax.plot(Vr, pr, lw=1, color='b',linestyle='dashed')[0]
line2 = ax.plot(Vr, pr, lw=1, color='b')[0]
line3 = ax.fill_between([],[],color='b',alpha=0.3)

areax1 = []
areay1 = []
line4 = ax.fill_between(areax1,areay1,color='orange',alpha=0.3)
areax2 = []
areay2 = []
line5 = ax.fill_between(areax2,areay2,color='orange',alpha=0.3)


textT = ax.text(0.7,0.74,'T$_r$ = {:.2f}'.format(Tr),transform=ax.transAxes)
textT1 = ax.text(0.7,0.68,'T = {:.0f}'.format(Tr*Tc),transform=ax.transAxes)

textflg = ax.text(1,1,'')
textfl = ax.text(1,1,'')
textg = ax.text(1,1,'')

# info box
infobox = ''
infobox += 'gas: '+gas + '\n' 
infobox += 'p$_c$: ' + "{:.2f}".format(pc/1e5) + ' (bar)\n' 
infobox += 'T$_c$: ' + "{:.0f}".format(Tc) + ' (K)\n' 
infobox += 'V$_c$: ' + "{:.0f}".format(Vc*1e6) + ' (cm$^3$/mol)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.7,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)



Tup = 350
Tdown = 257
Tsteps = 1200
new = True

def animate(k):
   global line3,line4,line5
   global line1,line2
   global textflg
   global textfl
   global textg
   global new
   T = Tup - (Tup-Tdown)*k/Tsteps 
   Tr = T/Tc
   textT.set_text('T$_r$ = {:.2f}'.format(Tr))
   textT1.set_text('T = {:.0f}'.format(T)+' K')
   
   if Tr<1:
     line1.set_ydata(vdw(Tr,Vr))
     line2.set_ydata(vdw_maxwell(Tr,Vr)[0])
     
     if new:
         areax1.append(1)
         areay1.append(1)
         areax2.append(1)
         areay2.append(1)
         new = False
     Vmin = vdw_maxwell(Tr, Vr)[1]
     pmin = vdw(Tr,Vmin)
     areax1.append(Vmin)
     areay1.append(pmin)
     line4.remove()
     line4 = ax.fill_between(areax1,areay1,color='orange',alpha=0.3,edgecolor=None)  
     
     Vmax = vdw_maxwell(Tr, Vr)[2]
     pmax = vdw(Tr,Vmax)
     areax2.append(Vmax)
     areay2.append(pmax)
     line5.remove()
     line5 = ax.fill_between(areax2,areay2,color='orange',alpha=0.3,edgecolor=None)
     
     line3.remove()
     line3 = ax.fill_between(Vr,vdw(Tr,Vr),vdw_maxwell(Tr, Vr)[0],color='blue',alpha=0.3)

     textflg.remove()
     textflg = ax.text(0.5*(Vmin+Vmax),pmin-0.05,'liquid/gas',horizontalalignment='center',verticalalignment='center',fontsize=8)
 
     textfl.remove()
     textfl = ax.text(0.5*(Vmin),pmin,'liquid',horizontalalignment='center',verticalalignment='center',fontsize=8)
 
     textg.remove()
     textg = ax.text(0.5*(Vmax+4),pmin,'gas',horizontalalignment='center',verticalalignment='center',fontsize=8)
 
    
   else:
     line1.set_ydata(prz)
     line2.set_ydata(vdw(Tr,Vr))
     
   
anim = animation.FuncAnimation(fig,animate,interval=1,frames=Tsteps)
anim.save('vanderWaals.mp4',fps=25,dpi=300)