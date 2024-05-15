#
# Icebow
# Examples for physics lectures
# Achim von Keudell
# Code from Peter M. Attia following a Fortran example in the appendix of 
# Bratt & Faulkner "Electrochemical methods fundamentals and applicatons" Wiley, 2001 (in Zotero)
# converted to python
# Ruhr University Bochum, 2024
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# INDEPENDENT VARIABLES
C = 1.0    # [=] mol/L, initial concentration of O. Default = 1.0
D = 1E-5   # [=] cm^2/s, O & R diffusion coefficient. Default = 1E-5
etai = +0.3   # [=] V, initial overpotential (relative to redox potential). Default = +0.2
etaf = -0.3   # [=] V, final overpotential (relative to redox potential). Default = -0.2
v = 2E-3   # [=] V/s, sweep rate. Default = 1E-3
n = 1.0    # [=] number of electrons transferred. Default = 1
alpha = 0.5   # [=] dimensionless charge-transfer coefficient. Default = 0.5
k0 = 1E-2   # [=] cm/s, electrochemical rate constant. Default = 1E-2
kc = 1E-3   # [=] 1/s, chemical rate constant. Default = 1E-3
T = 298.15  # [=] K, temperature. Default = 298.15

# PHYSICAL CONSTANTS
F = 96485   # [=] C/mol, Faraday's constant
R = 8.3145  # [=] J/mol-K, ideal gas constant
f = F / (R * T)  # [=] 1/V, normalized Faraday's constant at room temperature

# SIMULATION VARIABLES
L = 800    # [=] number of iterations per t_k (pg 790). Default = 500
DM = 0.45  # [=] model diffusion coefficient (pg 788). Default = 0.45

# DERIVED CONSTANTS
tk = 2 * (etai - etaf) / v    # [=] s, characteristic exp. time (pg 790). In this case, total time of fwd and rev scans
Dt = tk / L               # [=] s, delta time (Eqn B.1.10, pg 790)
Dx = np.sqrt(D * Dt / DM)      # [=] cm, delta x (Eqn B.1.13, pg 791)
j = int(np.ceil(4.2 * L**0.5) + 5)  # number of boxes (pg 792-793). If L~200, j=65

# REVERSIBILITY PARAMETERS
ktk = kc * tk               # dimensionless kinetic parameter (Eqn B.3.7, pg 797)
km = ktk / L               # normalized dimensionless kinetic parameter (see bottom of pg 797)
Lambda = k0 / (D * f * v)**0.5     # dimensionless reversibility parameter (Eqn 6.4.4, pg. 236-239)

# CHEMICAL REVERSIBILITY WARNING
if km > 0.1:
    print('Warning: k_c*t_k/l equals', km,
          ', which exceeds the upper limit of 0.1 (see B&F, pg 797)')

# PRE-INITIALIZATION
C /= 1000           # Convert C from mol/L to mol/cm3
k = np.arange(L+1)                # time index vector
t = Dt * k             # time vector
eta1 = etai - v * t      # overpotential vector, negative scan
eta2 = etaf + v * t      # overpotential vector, positive scan
eta = np.concatenate((eta1[eta1 > etaf], eta2[eta2 <= etai]))  # overpotential scan, both directions
Enorm = eta * f          # normalized overpotential
kf = k0 * np.exp(-alpha * n * Enorm)  # [=] cm/s, fwd rate constant (pg 799)
kb = k0 * np.exp((1 - alpha) * n * Enorm)  # [=] cm/s, rev rate constant (pg 799)

O = np.ones((L+1, j)) * C  # [=] mol/cm^3, concentration of O
R = np.zeros((L+1, j))      # [=] mol/cm^3, concentration of R
JO = np.zeros(L+1)          # [=] mol/cm^2-s, flux of O at the surface

#km = 0

# START SIMULATION
# i1 = time index. i2 = distance index
for i1 in range(L):
    # Update bulk concentrations of O and R
    # simple diffusion code
    for i2 in range(1, j-1):
        O[i1+1, i2] = O[i1, i2] + DM * (O[i1, i2+1] + O[i1, i2-1] - 2 * O[i1, i2])

        R[i1+1, i2] = R[i1, i2] + DM * (R[i1, i2+1] + R[i1, i2-1] - 2 * R[i1, i2]) \
            - km * R[i1, i2]

    # Update flux
    # calculate surface flux
    JO[i1+1] = (kf[i1+1] * O[i1+1, 2] - kb[i1+1] * R[i1+1, 2]) / (1 + Dx/D * (kf[i1+1] + kb[i1+1]))

    # Update surface concentrations
    O[i1+1, 0] = O[i1+1, 1] - JO[i1+1] * (Dx / D)
    R[i1+1, 0] = R[i1+1, 1] + JO[i1+1] * (Dx / D) - km * R[i1+1, 0]

# Calculate current density, Z, from flux of O
Z = -n * F * JO * 1000  # [=] A/cm^2 -> mA/cm^2, current density

# PLOT RESULTS
# Sometimes length(eta) = length(Z) + 1. If this is the case, truncate last value
if len(eta) > len(Z):
    eta = eta[:-1]

fig, ax = plt.subplots(2,1,figsize=(5,5))  
fig.tight_layout()
ax[0].plot(eta, Z)
lineZ = ax[0].plot(eta[0],Z[0],marker='o',color='red')[0]
ax[0].set_xlabel('overpotential (V)')
ax[0].set_ylabel('current density (mA/cm^2)')

infobox = ''
infobox += 'Sweep rate v: '+ "{:.3f}".format(v) + 'V/s'
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[0].text(0.02,0.95,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax[0].transAxes) 

lineR = ax[1].plot(R[100,:],label='R')[0]
lineO = ax[1].plot(O[100,:],label='O')[0]
ax[1].legend()
ax[1].fill_between([-5,0],[0.0015,0.0015],y2=0,fc='lightgrey')
ax[1].text(0.01,0.5,'electrode', fontsize=8,rotation=90,verticalalignment='top',transform=ax[1].transAxes)   
ax[1].set_xlabel('z (a.u.)')
ax[1].set_ylabel('concentrations O,R (a.u.)')
ax[1].set_yticks([])
ax[1].set_ylim(0,0.0015)
ax[1].set_xlim(-5,100)


infobox = ''
infobox += 'Reaction at electrode: O + e$^-$ > R'#+ "{:.0f}".format(NBreflectionsinside)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[1].text(0.1,0.95,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax[1].transAxes)   

def animate(k):
     
     if k%50 == 0: print(k,' of ',L) 
     lineZ.set_xdata(eta[k])
     lineZ.set_ydata(Z[k])
     lineR.set_ydata(R[k,:])
     lineO.set_ydata(O[k,:])
     
        
anim = animation.FuncAnimation(fig,animate,interval=800,frames=L)
anim.save('CV.mp4',fps=25,dpi=300)
