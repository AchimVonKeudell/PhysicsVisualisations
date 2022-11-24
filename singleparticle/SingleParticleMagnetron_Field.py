#
# Single particle movement in a 2'' magnetron
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Implementation of the Bfield by
# Julian Held
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.constants as const
from scipy.special import hyp2f1
import plotly.graph_objects as go
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot


el = -1.6e-19 # el charge
me = 9e-31 # el mass
mp = 1.6e-27 # proton mass

   

# --------------------------------------------------------------------------
# Define fields
# --------------------------------------------------------------------------
def Bfield(x):
    """ Returns Br, Bz in T for r, z position (in m).
    Reconstructed for the Thin Films consulting 2" UHV magnetrons.
    Reconstruction performed by Dennis Krueger (TET).
    z = 0 is the target surface, if traget_thickness is set correctly.
    z = 0 is the magnetron surface if target_thickness=0

    If used, cite: Krueger et al, Phys. Plasmas 25, 061207 (2018)"""
    target_thickness=3e-3
    R = 20.57870655967374e-3 #radius of the ring magnet (m)
    dC = 16.68888450585952e-3 #location of center magnet below the target (m)
    dR = 8.240298921866946e-3 #location of the ring magnet below the target (m)
    mC = -6.07855 #magnetic dipole moment of the center magnet (Am^2)
    mR = 10.2386 #magnetic dipole moment of the ring magnet (Am^2)

    r = np.sqrt(x[0]**2+x[1]**2)
    z = x[2]
    z = z + target_thickness

    Br = (3. * const.mu_0 * mC/(4.*np.pi) * (r*(z + dC)) / (((z+dC)**2 + r**2)**(5./2.))
    + 3. * const.mu_0 * mR / (4.*np.pi)  * (r*(z + dR))
    / (((z+dR)**2 + r**2 +R**2)**(5./2.))
    * hyp2f1(3./4., 5./4., 1., (4. * r**2 * R**2)
    / (((z+dR)**2 + r**2 + R**2)**2))
    - 15. * const.mu_0 * mR / (8.*np.pi) *  (r*R**2*(z+dR) *( (z+dR)**2-r**2+R**2)
    /(((z+dR)**2 + r**2 + R**2))**(9./2.))
    * hyp2f1(7./4., 9./4., 2., (4. * r**2 * R**2)
    / (((z+dR)**2+r**2+R**2)**2)))


    Bz = (const.mu_0 * mC/(4.*np.pi) * (2.*(z+dC)**2 - r**2)/((z+dC)**2 + r**2)**(5./2.)
       + const.mu_0 * mR/(4.*np.pi) * (2.*(z+dR)**2 - r**2 - R**2)
       / ((z+dR)**2 + r**2 + R**2)**(5./2.)
       * hyp2f1(3./4., 5./4., 1., (4. * r**2 * R**2)
       / ( ((z+dR)**2 + r**2 + R**2)**2) )
       + 15. * const.mu_0 * mR / (4.*np.pi) * (r**2*R**2*(z+dR)**2)
       / (((z+dR)**2 + r**2 + R**2)**(9./2.))
       * hyp2f1(7./4., 9./4., 2., (4. * r**2 * R**2)/((z+dR)**2+r**2+R**2)**2))

    if x[0] == 0:
       if x[1]>=0:
          phi = np.pi/2
       else:
          phi = 3/2*np.pi
    else:
       if x[0] >=0 and x[1]>=0: 
         phi = np.arctan(abs(x[1]/x[0]))
       elif x[0] <=0 and x[1]>=0: 
         phi = np.pi-np.arctan(abs(x[1]/x[0]))           
       elif x[0] <=0 and x[1]<=0: 
         phi = np.pi+np.arctan(abs(x[1]/x[0]))           
       elif x[0]>=0 and x[1]<=0: 
         phi = 2*np.pi-np.arctan(abs(x[1]/x[0]))           
         
    
    Bx = Br*np.cos(phi)
    By = Br*np.sin(phi)        
        
    return np.array([Bx,By,Bz]) # return field due to magnetic bottle


# --------------------------------------------------------------------------
# Setup Plot
# --------------------------------------------------------------------------
sa = 0.2/1e-3 # plot axis scale in mm
hvsw = 1
gxy = 10
gz = 10


# define field arrows in the plot in mm
# quiver can only generate arrows with different length by
# using individual axp.quiver commands
xgrid, ygrid, zgrid = np.meshgrid(np.arange(-sa,sa+sa/(gxy),2*sa/(gxy)),
                                  np.arange(-sa,sa+sa/(gxy),2*sa/(gxy)),
                                  np.arange(0,hvsw*2*sa+hvsw*sa/gz,hvsw*2*sa/gz))

ubfeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
vbfeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
wbfeld = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
lengthb = np.zeros((len(xgrid[:,0,0]),len(xgrid[0,:,0]),len(xgrid[0,0,:])))
for i in range(len(xgrid[:,0,0])):
    for j in range(len(xgrid[0,:,0])):
        for m in range(len(xgrid[0,0,:])):
            # real points in m, xgrid points in mm
            xb = xgrid[i,j,m]*1e-3
            yb = ygrid[i,j,m]*1e-3
            zb = zgrid[i,j,m]*1e-3
            ubfeld[i,j,m], vbfeld[i,j,m], wbfeld[i,j,m] = Bfield([xb,yb,zb])
            lengthb[i,j,m] = np.sqrt( Bfield([xb,yb,zb])[0]**2 + 
                                      Bfield([xb,yb,zb])[1]**2 + Bfield([xb,yb,zb])[2]**2)


fig = go.Figure(data=go.Streamtube(x=xgrid,y=ygrid,z=zgrid,u=ubfeld,v=vbfeld,w=wbfeld,
                                   starts=dict(x=xgrid[:,0,0],y=ygrid[0,:,0],z=zgrid[:,:,0])))
#axp = plt.axes(projection='3d')
#axp.set(xlabel="x (mm)",ylabel="y (mm)",zlabel="z (mm)")
#axp.set(xlim=(-sa,sa),ylim=(-sa,sa),zlim=(0,2*hvsw*sa))
plot(fig)


 