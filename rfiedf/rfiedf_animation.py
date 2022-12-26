#
# IEDF in einer RF Randschicht
# Examples for plasma pyhsics lectures
# original code by Simon Kreuznacht
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

model = 'nocollisions'
model = 'collisions'

if model == 'nocollisions':
   animationfile = 'rfiedfnocollisions.mp4'
   num = 0 # 1e6 # collision frequency
   #Schichtmodel='Matrix' # Randschichtmodell
   Schichtmodel='Child-Langmuir' # Randschichtmodell
   mion = 4 # in amu
if model == 'collisions':
   animationfile = 'rfiedfcollisions.mp4'
   num = 4e6 # 1e6 # collision frequency
   #Schichtmodel='Matrix' # Randschichtmodell
   Schichtmodel='Child-Langmuir' # Randschichtmodell
   mion = 4 # in amu

# ----------------------------------------------------------------------------
#
# Monte Carlo Simulation der IEDF an der Elektrode in einer RF Randschicht 
#
# ----------------------------------------------------------------------------

el = 1.602176634e-19 # el charge
amu = 1.6605390666e-27 # proton mass
kB = 1.380649e-23 # Boltzmann const.

Te = 3 # electron temperature
smax = 0.008 # in m
V0 = 700 # Peak to peak
fRF = 13.6 # in MHz
omegaRF = 2*np.pi*fRF*1e6 # angular frequency

dt = 1e-10 # Schrittweite der Zeitschritte
NBParticles = 500
qion = 1 # in Elementarladungen
Ewidth = 5 # Width energy bin in eV
Emax = 500 # Max energy scale
nup = dt*num # probaility of collision per time step

# store IEDF
IEDF = np.zeros(int(Emax/Ewidth))
IEDFNorm = np.zeros(int(Emax/Ewidth))
Escan = np.linspace(0,Emax,int(Emax/Ewidth))

# setup plot
fig, axp = plt.subplots(1,1,figsize=(8,4.5))
line1 = axp.plot(Escan,IEDFNorm,color='b',lw=1, label ='')[0]
axp.set(xlabel="energy (eV)",ylabel="N(eV)/N_max")
axp.set(xlim=(0,Emax),ylim=(0,1))

# setup info box
infobox = ''
infobox += 'Model: ' + Schichtmodel + ' \n' 
infobox += 'V0: ' + "{:.0f}".format(V0) + ' (V)\n' 
infobox += 'fRF: ' + "{:.2f}".format(fRF) + ' (MHz)\n' 
infobox += 'Te: ' + "{:.1f}".format(Te) + ' (eV)\n' 
infobox += 'mion: ' + "{:.0f}".format(mion) + ' (amu)\n' 
infobox += 'num: ' + "{:.2e}".format(num) + ' (1/s)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
axp.text(0.05*Emax,0.95,infobox, fontsize=8,bbox=props,verticalalignment='top')

# Initialize Simulation
Particles = 0

def animate(k):
    
    global IEDF,Escan,Particles,IEDFNorm
    
    print(k,' of ',NBParticles)
    microsteps = 100 # calculate 100 particles before image is generated
    for k in range(microsteps):
        Particles += 1
        PhiStart = np.random.rand()*2*np.pi # phase of entry in the sheath
        x = smax # ion start at maximum sheath extension    
        # ion start with Bohm velocity
        v =-np.sqrt(Te*el/(mion*amu)) 
        t = 0 # start at t=0
        a = 0 # start without Force
        EFeld=0
        while x>0:
            t += dt
            x += v*dt + 1/2*dt**2*a # new position        
            # calculation of sheath thickness and acceleration
            if Schichtmodel == 'Matrix':
                # sheath thickness between 0 and smax
                s = smax/2*(1-np.sin(omegaRF*t+PhiStart)) 
                if x<=s and s>=0:     
                    # when the sheath collapsed, the volate is zero
                    # maximal extension corresponds to peak-to-peak voltage
                    VSchicht = V0/4*(1-np.sin(omegaRF*t + PhiStart))**2 
                    EFeld = 2*VSchicht/(s**2)*(x-s) # Matrix Model
                else:
                    EFeld=0                
            elif Schichtmodel == 'Child-Langmuir':
                s = smax/2*(1-np.sin(omegaRF*t + PhiStart));
                if x<=s and s>0:                
                    # when the sheath collapsed, the volate is zero
                    # maximal extension corresponds to peak-to-peak voltage
                    VSchicht = V0/4*(1-np.sin(omegaRF*t + PhiStart))**2 
                    EFeld = -4/3*VSchicht/(s**(4/3))*(abs(x-s))**(1/3) # Child-Langmuir Model
                else:
                    EFeld=0
            a = qion*el*EFeld/(mion*amu) # acceleration
            v +=a*dt # propagate velocity
        
            # CX collisions
            if np.random.rand()<=nup:
                v=0 
        
        Eion=1/2*mion*amu*v**2/el # in eV
        Ebin=round(Eion/Ewidth) # select bin
        if Ebin>=1 and Ebin<=Emax/Ewidth:
            IEDF[Ebin] += 1 # increase number of ions in Bin
   
    IEDFNorm = IEDF/max(IEDF) # normalize
    
    line1.set_xdata(Escan)
    line1.set_ydata(IEDFNorm) 
    fig.suptitle('Particles: ' + "{:.0f}".format(Particles))
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=NBParticles)
anim.save(animationfile,fps=25,dpi=300)

plt.plot(Escan,IEDFNorm)
plt.xlabel('Ion energy [eV]')
plt.ylabel('N(E) [normalized]')
