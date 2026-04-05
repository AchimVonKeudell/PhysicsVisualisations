# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 06:23:32 2023

@author: Achim
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation
from scipy import optimize

kb = 1.38e-23
el = 1.6e-19
h = 6.6e-34
m = 9.31e-31
ener = np.linspace(1e-9,10,1000)

model = 'EF Temp variation'
model = 'Density variation'

if model == 'EF Temp variation':
    ef = 4
    mu = 0
    steps = 1000
    T = 300
    animationname = 'qstatistics_temp.mp4'
    
if model == 'Density variation':
    ef = 4
    ef0 = 4
    mu = 0
    steps = 100
    T = 1000
    animationname = 'qstatistics_density.mp4'

def de(ener):
    return 8*np.pi/(h**3)*np.sqrt(2*m**3)*np.sqrt(ener*el)*1e-27*el

def fd(en,ef,T):
    return de(en)*1/(np.exp((en-ef)*el/(kb*T))+1)

def mb(en,ef,T):
    return de(en)*np.exp(-(en)*el/(kb*T))

def be(en,mu,T):
    #z = np.exp(mu/(kb*T))
    return de(en)*1/(np.exp((en-mu)*el/(kb*T))-1)#+z/(1-z)


# ----------------------------------------------------------------------------
# Initialize plots
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(1,1,figsize=(8,4.5))

if model == 'EF Temp variation':
  line1 = ax.plot(ener,fd(ener,ef,T),label='F.D.')[0]
  line3 = ax.plot([ef,ef],[0,20],label='E_F',color='red',linestyle='dashed')[0]

  ax.set(xlabel='energy (eV)',ylabel='delta n(E) per nm^3 and eV',xlim=(0,10),ylim=(0,20))
  ax.legend()

if model == 'Density variation':
  line1 = ax.plot(ener,fd(ener,ef,T),label='F.D.')[0]
  line2 = ax.plot(ener,be(ener,mu,T),label='B.E.')[0]
  line3 = ax.plot(ener,mb(ener,ef,T),label='M.B.')[0]

  ax.set(xlabel='energy (eV)',ylabel='delta n(E) per nm^3 and eV',xlim=(0,10),ylim=(0,20))
  ax.legend()

normfd0 = integrate.quad(fd,1e-5,20,args=(ef,1e-9))[0]
#print('FDnorm: ',norm0)
normbe0 = integrate.quad(be,1e-5,20,args=(mu,1e-9))[0]
#print('BEnorm: ',normbe0)


#print(norm)

#norm = integrate.quad(be,0,20,args=(mu,T))
#print(norm)

# --------------------------------------------------------------------------
# Main loop
# --------------------------------------------------------------------------
def animate(k):
  global pos,p,t,time,ef,mu,normfd0,normbe0

  if model == 'EF Temp variation':
    Tloop = T + k/1000*(5000-T)
    print('Temp: '+ "{:.2f}".format(Tloop)+' K ')
  
    # F.D distribution
    # find ef distribution so that the integral is constant.
    efermi = ef
    def f(eft):
       return (integrate.quad(fd,0,20,args=(eft,Tloop))[0] - normfd0)**2    
    sol = optimize.minimize(f,efermi, method='Nelder-Mead',tol=1e-8)   
    #print('set EF to:')
    #print(sol.x)
    ef = sol.x[0]
    #print(integrate.quad(fd,0,20,args=(ef,Tloop))[0])
    line1.set_ydata(fd(ener,ef,Tloop)) 
    line3.set_xdata([ef,ef]) 
  
    fig.suptitle('Temperature: ' + "{:.0f}".format(Tloop)+' K\n'+'E_F: '+"{:.3f}".format(ef)+' eV')
   
  if model == 'Density variation':

    efermi = ef0 - k/100*ef0+0.01
    Tloop = T + k/100*(5000-T)
    normfd0 = integrate.quad(fd,1e-5,20,args=(efermi,1e-9))[0]

    print('E_F: '+ "{:.2f}".format(efermi)+' eV ',k)

    def f(eft):
       return (integrate.quad(fd,0,20,args=(eft,Tloop))[0] - normfd0)**2    
    sol = optimize.minimize(f,efermi, method='Nelder-Mead',tol=1e-8)   
    #print('set EF to:')
    #print(sol.x)
    ef = sol.x[0]
    if ef<0: ef = 0
    #print(integrate.quad(fd,0,20,args=(ef,Tloop))[0])
    line1.set_ydata(fd(ener,ef,Tloop)) 
   
    # B.E distribution
    # find ef distribution so that the integral is constant.
    #chempot = mu
    #def f(chempo):
    #  return (integrate.quad(be,1e-9,20,args=(chempo,Tloop))[0] - normfd0)**2    
    #sol = optimize.minimize(f,chempot, method='Nelder-Mead',tol=1e-8)   
    #print(sol.x)#
    #mu = sol.x[0]
    #print('set mu to: ',mu)
    #print(integrate.quad(be,0,20,args=(mu,Tloop))[0])
    
    line2.set_ydata(be(ener,0,Tloop)) 
    line3.set_ydata(mb(ener,0,Tloop)) 
    #line2.set_ydata(de(ener)) 
  
      
    fig.suptitle('Temperature: ' + "{:.0f}".format(Tloop)+' K\n'+'E_F: '+"{:.3f}".format(ef)+' eV,'+'mu: '+"{:.3f}".format(mu)+' eV')
    
    
anim = animation.FuncAnimation(fig,animate,interval=1,frames=steps)
anim.save(animationname,fps=25,dpi=300)
    