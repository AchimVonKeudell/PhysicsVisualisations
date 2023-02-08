#
# Beam bending
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
# basis was python code from Justin Verdirame
# https://mechanicsandmachines.com/?p=705
#

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


#model = 'cantilever_forcecenter'
model = 'beamsupport_forcecenter'
#model = 'beamclamped_forcecenter' 
#model = 'beamclampedsliding_forcecenter' 
model = 'beamclamped_momentcenter'
#model = 'beamcsupport_momentcenter'
#model = 'beamclamped_weight'
#model = 'beamcsupport_weight'
#model = 'beamsupportasymmetric_forcecenter'




if model == 'cantilever_forcecenter': 
    bcs = 'cantilever'
    loadcase = 'pointforce'
    forceposition =  0.5
    maxforce = 10
    animationfile = 'bendingcantileverforcecenter.mp4'
    ymin = -0.03
    ymax = 0.005
if model == 'beamsupport_forcecenter': 
    bcs = 'simplysupported'
    loadcase = 'pointforce'
    forceposition =  0.5
    maxforce = 30
    animationfile = 'bendingbeamforcecenter.mp4'
    ymin = -0.03
    ymax = 0.005
if model == 'beamsupportasymmetric_forcecenter': 
    bcs = 'simplysupportedasymetric'
    loadcase = 'pointforce'
    forceposition =  0.5
    maxforce = 30
    animationfile = 'bendingbeamforcecenterasymmetric.mp4'
    ymin = -0.03
    ymax = 0.005
if model == 'beamclamped_forcecenter': 
    bcs = 'clamped-clamped'
    loadcase = 'pointforce'
    forceposition =  0.5
    maxforce = 100
    animationfile = 'bendingbeamclampedforcecenter.mp4'
    ymin = -0.03
    ymax = 0.005
if model == 'beamclampedsliding_forcecenter': 
    bcs = 'clamped-sliding'
    loadcase = 'pointforce'
    forceposition =  0.5
    maxforce = 15
    animationfile = 'bendingbeamclammpedslidingforcecenter.mp4'
    ymin = -0.03
    ymax = 0.005
if model == 'beamclamped_momentcenter': 
    bcs = 'clamped-clamped'
    loadcase = 'pointmoment'
    forceposition =  0.5
    maxforce = 30
    animationfile = 'bendingbeamclampedmomentcenter.mp4'
    ymin = -0.01
    ymax = 0.01
if model == 'beamcsupport_momentcenter': 
    bcs = 'simplysupported'
    loadcase = 'pointmoment'
    forceposition =  0.5
    maxforce = 30
    animationfile = 'bendingbeamsupportmomentcenter.mp4'
    ymin = -0.01
    ymax = 0.01
if model == 'beamclamped_weight': 
    bcs = 'clamped-clamped'
    loadcase = 'uniformpres'
    forceposition =  0.5
    maxforce = 30
    animationfile = 'bendingbeamclampedweight.mp4'
    ymin = -0.002
    ymax = 0.002
if model == 'beamcsupport_weight': 
    bcs = 'simplysupported'
    loadcase = 'uniformpres'
    forceposition =  0.5
    maxforce = 30
    animationfile = 'bendingbeamsupportweight.mp4'
    ymin = -0.002
    ymax = 0.002
   
# Properties of a steel bar 10 mm x 5 mm    

E=200e9 #Young's modulus = 200 GPa
rho=7800. #density
L=1. #total length
b=0.010 #beam depth
h=0.005 #beam height
A=b*h #beam cross sectional area
I=b*h**3/12 #beam moment of inertia
N=40 #number of elements 
Nnodes=N+1 #number of nodes
ell=L/N #length of beam element
Nmodes=12 #number of modes to keep
mass=rho*L*A #total mass
Irot=1/12*mass*L**2 #rotary inertia
x=np.arange(0,ell+L,ell)
y = np.zeros(Nnodes)

K_e=E*I/ell**3*np.array([[12, 6*ell, -12, 6*ell],
                         [6*ell, 4*ell**2, -6*ell, 2*ell**2],
                         [-12, -6*ell, 12, -6*ell], #x2
                         [6*ell, 2*ell**2, -6*ell, 4*ell**2]]) #theta2


# Input boundary conditions here
#bcs = 'clamped-clamped'
#bcs = 'clamped-sliding'
#bcs = 'simplysupported'
#bcs = 'cantilever'

# Input load case here
#loadcase = 'uniformpres'
#loadcase = 'pointforce'
#loadcase = 'pointmoment'
#loadcase = 'unifromdistribmoment' 
#loadcase = 'uniformpres'
#loadcase = 'concmomenteachnode'

#forceposition =  0.5
#maxforce = 30

def setforcevector(force):
    global K,f,Nnodes
    
    f=np.zeros((2*Nnodes,1))

    
    if loadcase == 'pointforce':
        #print('Load: concentrated force at x=L/2')
        #load: concentrated force on center of beam
        #f[Nnodes-1]=force
        f[2*(int(forceposition*Nnodes))]=force
        return f
    elif loadcase == 'pointmoment':
        #print('Load: concentrated moment at x=L/2')
        #load: concentrated moment on center of beam
        f[Nnodes]=force 
        return f
    elif loadcase == 'uniformdistribmoment':
        #print('Load: uniformly distributed moment')
        #load: uniform distributed moment on entire beam
        m=1.0
        m_el=np.zeros((4,1))
        m_el[:,0] = np.array([-m,0,m,0])
        for i in range(1,N+1):
            f[2*(i-1)+1-1:2*(i-1)+4]=f[2*(i-1)+1-1:2*(i-1)+4]+m_el
        return f    
    elif loadcase == 'concmomenteachnode':
        #print('Load: Concentrated moment on each node')
        #load: concentrated moment on each node
        f[1:2*Nnodes:2]=force 
        return f
    elif loadcase == 'uniformpres':
        #print('Uniform pressure')
        #load: uniform distributed load on entire beam
        q=-rho*A*9.81*force/maxforce*(-1)
        q_el=np.zeros((4,1))
        q_el[:,0] = np.array([q*ell/2, q*ell**2/12, q*ell/2, -q*ell**2/12])
        for i in range(1,N+1):
            f[2*(i-1)+1-1:2*(i-1)+4]=f[2*(i-1)+1-1:2*(i-1)+4]+q_el
        return f    
    else:
        print('Unknown loading')

def boundarycondition():
    global K,f,Nnodes
    
    K=np.zeros((2*(N+1),2*(N+1))) 
    # M=zeros(2*(N+1),2*(N+1))
    for i in range(1,N+1):
        K[(2*(i-1)+1-1):(2*(i-1)+4), (2*(i-1)+1-1):(2*(i-1)+4)]=K[(2*(i-1)+1-1):(2*(i-1)+4)][:,(2*(i-1)+1-1):(2*(i-1)+4)] +K_e

    if bcs == 'clamped-clamped':
        #print('BCs: clamped-clamped')
        #clamped-clamped beam BCs
        K=np.delete(K,[2*Nnodes-2,2*Nnodes-1],0) #K[1:2,:]=[] 
        K=np.delete(K,[2*Nnodes-2,2*Nnodes-1],1) #K[:,1:2]=[] 
        f=np.delete(f,[2*Nnodes-2,2*Nnodes-1]) #f[1:2]=[] 
        K=np.delete(K,[0,1],0) #K[1:2,:]=[] 
        K=np.delete(K,[0,1],1) #K[:,1:2]=[] 
        f=np.delete(f,[0,1]) #f[1:2]=[]

    elif bcs == 'simplysupported':
        #print('BCs: simply supported')
        #simply supported beam BCs
        K=np.delete(K,[2*Nnodes-2],0) #K[1:2,:]=[] 
        K=np.delete(K,[2*Nnodes-2],1) #K[:,1:2]=[] 
        f=np.delete(f,[2*Nnodes-2]) #f[1:2]=[] 
        K=np.delete(K,[0],0) #K[1:2,:]=[] 
        K=np.delete(K,[0],1) #K[1:2,:]=[] 
        f=np.delete(f,[0]) #f[1:2]=[] 
   
    elif bcs == 'simplysupportedasymetric':
        #print('BCs: simply supported')
        #simply supported beam BCs
        K=np.delete(K,[2*Nnodes-2-20],0) #K[1:2,:]=[] 
        K=np.delete(K,[2*Nnodes-2-20],1) #K[:,1:2]=[] 
        f=np.delete(f,[2*Nnodes-2-20]) #f[1:2]=[] 
        K=np.delete(K,[0],0) #K[1:2,:]=[] 
        K=np.delete(K,[0],1) #K[1:2,:]=[] 
        f=np.delete(f,[0]) #f[1:2]=[]    
        
 
    elif bcs == 'cantilever':
        #print('BCs: cantilever')
        #cantilever beam BCs
        K=np.delete(K,[0,1],0) #K[1:2,:]=[] 
        K=np.delete(K,[0,1],1) #K[:,1:2]=[] 
        f=np.delete(f,[0,1]) #f[1:2]=[] 
 
    elif bcs == 'clamped-sliding':
        #print('BCs: clamped-sliding')
        #clamped-sliding beam BCs
        K=np.delete(K,[2*Nnodes-1],0) #K[1:2,:]=[] 
        K=np.delete(K,[2*Nnodes-1],1) #K[:,1:2]=[] 
        f=np.delete(f,[2*Nnodes-1]) #f[1:2]=[] 
        K=np.delete(K,[0,1],0) #K[1:2,:]=[] 
        K=np.delete(K,[0,1],1) #K[:,1:2]=[] 
        f=np.delete(f,[0,1]) #f[1:2]=[] 
    else:
        print('Unknown boundary conditions.')

# setup plot
fig, ax = plt.subplots(1,1,figsize=(8,4.5))
ax.set(ylabel='displacement [cm]',xlabel='position [m]',ylim=(ymin/1e-2,ymax/1e-2),xlim=(-0.1,1.1))
line = ax.plot(x,y)[0]

# info box
infobox = ''
infobox += 'loading: ' + loadcase +'\n'
infobox += 'boundary: ' + bcs + ''
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.1,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)



if loadcase == 'pointforce':
    linequiver = ax.quiver(x[int(forceposition*Nnodes)],0,0,-1,color='r',pivot='tip')
if loadcase == 'pointmoment':    
    xct = []
    yct = []
    for i in range(90):
        xc = x[int(Nnodes*forceposition)]+0.02*np.cos(i/90*2*np.pi)
        yc = 0.001/1e-2*np.sin(i/90*2*np.pi)
        xct.append(xc)
        yct.append(yc)        
    ax.plot(xct,yct,color='r',lw=0.5)    
    ax.quiver(x[int(forceposition*Nnodes)]+0.02,0,0,-1,color='r',pivot='middle',scale=30,width=0.01)
    ax.quiver(x[int(forceposition*Nnodes)]-0.02,0,0,1,color='r',pivot='middle',scale=30,width=0.01)
if loadcase == 'uniformpres':
    linequiver = ax.quiver(x,y,0,-1,color='r',pivot='tip',scale=20,width=0.005)
    
# display bc
if bcs == 'simplysupported':
    ax.plot(0,0,marker='o',markersize=5,color='g')
    ax.plot(x[Nnodes-1],0,marker='o',markersize=5,color='g')
if bcs == 'simplysupportedasymetric':
    ax.plot(0,0,marker='o',markersize=5,color='g')
    ax.plot(x[Nnodes-1-5],0,marker='o',markersize=5,color='g')
if bcs == 'clamped-clamped':
    ax.plot([-0.05,0],[0,0],lw=5,color='g')
    ax.plot([x[Nnodes-1],x[Nnodes-1]+0.05],[0,0],lw=5,color='g')
if bcs == 'cantilever':
    ax.plot([-0.05,0],[0,0],lw=5,color='g')
if bcs == 'clamped-sliding':
    ax.plot([-0.05,0],[0,0],lw=5,color='g')
    lineright = ax.plot([x[Nnodes-1],x[Nnodes-1]+0.05],[0,0],lw=1,linestyle='dashed',color='g')[0]


def animate(p):
    global K,f
    # increase and decrease  the force
    if p<100:
       force = -maxforce*p/200
    else:
       force = -maxforce*(200-p)/200
       
    f = setforcevector(force)
    boundarycondition()
    
    
    # solvinge system of equations
    dx_vec=np.linalg.solve(K,f)    
    
    if bcs == 'clamped-clamped':
        #print('BCs output: clamped-clamped')
        #clamped-clamped
        dx=np.hstack([0., dx_vec[0:2*Nnodes-5:2], 0.]) 
        dtheta=np.hstack([0., dx_vec[1:2*Nnodes-4:2], 0.])
    elif bcs == 'simplysupported':
        #print('BCs output: simply-supported')
        #simply-supported
        dx=np.hstack([0., dx_vec[1:2*Nnodes-4:2], 0.])
        dtheta=np.hstack([dx_vec[0:2*Nnodes-3:2], dx_vec[2*Nnodes-3]])
    elif bcs == 'simplysupportedasymetric':
        #print('BCs output: simply-supported')
        #simply-supported
        dx=np.hstack([0., dx_vec[1:2*Nnodes-4:2], 0.])
        dtheta=np.hstack([dx_vec[0:2*Nnodes-3:2], dx_vec[2*Nnodes-3]])
    elif bcs == 'cantilever':
        #print('BCs output: cantilever')
        #cantilever
        dx=np.hstack([0., dx_vec[0:2*Nnodes-2:2]])
        dtheta=np.hstack([0., dx_vec[1:2*Nnodes-1:2]])
    elif bcs == 'clamped-sliding':
        #print('BCs output: clamped-sliding')
        #clamped-sliding beam BCs
        dx=np.hstack([0., dx_vec[0:2*Nnodes-3:2]])
        dtheta=np.hstack([0., dx_vec[1:2*Nnodes-4:2], 0.])
    else:
        print('Output: Unknown boundary conditions.')

    #print(dx)
    #print(x)
    #line.set_xdata(x)
    line.set_ydata(dx/1e-2)
    if loadcase =='pointforce':
        linequiver.set_offsets([x[int(forceposition*Nnodes)],dx[int(forceposition*Nnodes)]/1e-2])
        fig.suptitle('Force: ' + "{:.1f}".format(force)+' N')
    if loadcase == 'uniformpres':        
        linequiver.set_offsets(np.c_[x,dx/1e-2])
        fig.suptitle('Force: weight')
    if loadcase == 'pointmoment':
        fig.suptitle('Torque: ' + "{:.1f}".format(force)+' Nm')
    if bcs == 'clamped-sliding':
        lineright.set_ydata([dx[Nnodes-1]/1e-2,dx[Nnodes-1]/1e-2])
    
        
    
    
    

anim = animation.FuncAnimation(fig,animate,interval=1,frames=200)
anim.save(animationfile,fps=25,dpi=300)