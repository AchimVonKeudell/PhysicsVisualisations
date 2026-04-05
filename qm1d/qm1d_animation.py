#
# Wave packet 1d solution
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#
# Initial Code von Louis Sharma, Paris
# https://github.com/louishrm/Quantum-Tunneling/blob/main/QM%20Tunnelling.ipynb
#

import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation

model = 'Barrier'
model = '2Barrier'
model = 'Coulomb'

if model == 'Barrier':
  filename = 'barrier.gif'    
  N_Grid = 750
  L = 750
  a = 320    
  V0 = 0.1
  w = 40
  x0 = 100
  k0 = 0.4
  sigma = 15
  tend = 1500
if model == '2Barrier':
  filename = '2barrier.gif'    
  N_Grid = 750
  L = 750
  a = 320 
  a2 = 450   
  V0 = 0.1
  w = 50
  w2 = 30
  x0 = 100
  k0 = 0.4
  sigma = 15
  tend = 1500
if model == 'Coulomb':
  filename = 'coulomb.gif'    
  N_Grid = 750
  L = 750
  a = 320    
  V0 = 0.4
  w = 300
  x0 = 100
  k0 = 0.4
  sigma = 15
  tend = 1500
  
  

t = np.linspace(0.,tend,1000)      
x = np.linspace(0,L,N_Grid+1) #grid of points
dx = x[1]-x[0] #grid point spacing or 'discrete' analogue of the differential length
        
              
def integral(f,axis = 0):
            """This function allows us to approximate integrals in discrete space"""
            return np.sum(f*dx, axis = axis)
        
        
Psi0  = np.exp( -1/2* (x[1:-1]-x0)**2/sigma**2) *np.exp(1j*k0*x[1:-1]) 
#use this range for x because as mentionned, we need the wavefunction to be 0 at the endpoints of the grid. 
        
        
#normalise the initial state
norm  = integral(np.abs(Psi0)**2)
Psi0 = Psi0/np.sqrt(norm)
        
#kinetic energy
T = -1/2 * 1/dx**2 * (np.diag(-2*np.ones(N_Grid-1))+ np.diag(np.ones(N_Grid-2),1)+ np.diag(np.ones(N_Grid-2),-1))
        
#potential as a flat array
if model == 'Barrier':
  V_flat = np.array([V0 if a< pos < a+w else 0 for pos in x[1:-1]])
if model == 'Coulomb':
  V_flat = np.array([V0/(pos-a+0.01) if a< pos < a+w else 0 for pos in x[1:-1]])
if model == '2Barrier':
  V_flat = np.array([V0 if a< pos < a+w or a2< pos < a2+w2 else 0 for pos in x[1:-1]])
        
#potential energy as a diagonal matrix
V = np.diag(V_flat)
        
#Hamiltonian
H = T+V
        
fig = plt.figure(figsize = (10,6))
ax = plt.axes(xlim=(0, L), ylim=(-0.25, 0.25))
line, = ax.plot([], [], lw=2)
ax.fill_between(x[1:-1],V_flat, label = '$V(x)$',facecolor='lightblue')
#ax.set_title('Gaussian wave packet with a potential barrier', fontsize = 20)
line1, = ax.plot(x[1:-1],np.zeros(N_Grid-1),lw=2,color="orange", label = '$\Re(\psi)$')
line2, = ax.plot(x[1:-1],np.zeros(N_Grid-1),lw=2,color="blue", label = '$\Im(\psi)$')
ax.legend(fontsize = 15)
#ax.set_xlabel('$x$', fontsize = 15)
ax.axis('off')
line1.set_data([],[])  
line2.set_data([], [])
          

def integral(f,axis = 0):
            """This function allows us to approximate integrals in discrete space"""
            return np.sum(f*dx, axis = axis)
            
            
#get eigenvalues and eigenvectors and normalise
E, psi = np.linalg.eigh(H)
psi = psi.T
norm = integral(np.abs(psi)**2)
psi = psi/np.sqrt(norm)
Etotal = 0

#get expansion coeffs
c_n = np.zeros_like(psi[0], dtype=complex)
for j in range(0, N_Grid-1):

            c_n[j] = integral(np.conj(psi[j]) * Psi0) #for each eigenvector, compute the inner product
            Etotal += c_n[j]*E[j]
            
print(Etotal)            
        
#solve the eigenvalue problem and get the time-dependent wavefunction   
def animation(k):
        
        print(k)
        tt = t[k]
        
   
        #get a function that returns the time dependent wavefunction
        def Psi(t):
            
            return psi.T@(c_n*np.exp(-1j*E*t))


        y1 = np.real(Psi(tt))
        y2 = np.imag(Psi(tt))
        line1.set_data(x[1:-1],y1)  
        line2.set_data(x[1:-1], y2)


            
ani = FuncAnimation(fig, animation, frames=len(t), interval=20, blit=False)
ani.save(filename)
                    
#wavepacket = Gaussian_Wave(750,750, 420,0.1,50,100,0.4,15, np.linspace(0.,2500,1000))
#Psi = wavepacket.animation()            
            
            
            
            
        
        
        