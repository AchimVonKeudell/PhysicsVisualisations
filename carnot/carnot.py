# 
# Carnotn
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2024
#
# parts of the code from Advances in Engineering Software Volume 172, October 2022, 103186
#
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
import matplotlib.cm as cm


NBprocess = 200

model = 'Stirling'
#model = 'Carnot'
#model = 'Otto'
#model = 'Diesel'

if model == 'Stirling':
    animationfile = 'stirling.gif'
    process1 = 'adiabatic'
    process2 = 'isochor'
    process3 = 'adiabatic'
    process4 = 'isochor'
    NBsteps = 4*NBprocess
    Tscale = 2000
    kappa = 1.4
    compression = 5
if model == 'Carnot':
    animationfile = 'carnot.gif'
    process1 = 'isotherm'
    process2 = 'adiabatic'
    process3 = 'isotherm'
    process4 = 'adiabatic'
    NBsteps = 4*NBprocess
    Tscale = 800
    kappa = 1.4
    compression = 5
if model == 'Otto':
    animationfile = 'otto.gif'
    process1 = 'adiabatic'
    process2 = 'isochor'
    process3 = 'adiabatic'
    process4 = 'isochor'
    NBsteps = 4*NBprocess
    Tscale = 2000
    kappa = 1.4
    compression = 5
if model == 'Diesel':
    animationfile = 'diesel.gif'
    process1 = 'adiabatic'
    process2 = 'isobar'
    process3 = 'adiabatic'
    process4 = 'isochor'
    NBsteps = 4*NBprocess
    Tscale = 5000
    kappa = 1.4
    compression = 5
#~~~~~~~~~~~~~~~~~~~
#   Carnot Cycle
#~~~~~~~~~~~~~~~~~~~
def carnot(p_min,p_max,v_max,r,gma):
    '''This function prints carnot cycle
    arguments are as follows:
    p_min: minimum pressure
    p_max: Maximum pressure
    v_max: Maximum volume
    r: compression ratio
    gma: Adiabatic exponent
    The order of arguments is:
    p_min,p_max,v_max,r,gma
    
    NOTE: After calling the function you have to type:
    show()
    '''
    p1=p_min
    v1=v_max
    v2=v1/r
    c1=p1*v1
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 1-2
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v=linspace(v1,v2,NBprocess)
    p=c1/v
    step1p = p
    step1v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 2-3
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p3=p_max
    c2=c1*v2**(gma-1.)
    v3=(c2/p3)**(1/gma)
    v=linspace(v2,v3,NBprocess)
    p=c2/v**gma
    step2p = p
    step2v = v
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 3-4
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    c3=p3*v3
    c4=p1*v1**gma
    v4=(c4/c3)**(1/(gma-1.))
    v=linspace(v3,v4,NBprocess)
    p=c3/v
    step3p = p
    step3v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 3-4
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v=linspace(v4,v1,NBprocess)
    p=c4/v**gma
    step4p = p
    step4v = v

    return step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v
    
#~~~~~~~~~~~~~~~~~~~
#    Otto Cycle
#~~~~~~~~~~~~~~~~~~~

def otto(p_min,p_max,v_max,r,gma):
    '''This function prints Otto cycle
    arguments are as follows:
    p_min: minimum pressure
    p_max: Maximum pressure
    v_max: Maximum volume
    r: compression ratio
    gma: Adiabatic exponent
    The order of arguments is:
    p_min,p_max,v_max,r,gma
    
    NOTE: After calling the function you have to type:
    show()
    '''
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 1-2
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p1=p_min
    v1=v_max
    v2=v1/r
    c1=p1*v1**gma
    v=linspace(v1,v2,NBprocess)
    p=c1/v**gma
    step1p = p
    step1v = v
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 2-3
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p3=p_max
    v3=v2
    p2=c1/v2**gma
    p=linspace(p2,p3,NBprocess)
    v=NBprocess*[v3]
    step2p = p
    step2v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 3-4
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    c2=p3*v3**gma
    v4=v1
    v=linspace(v3,v4,NBprocess)
    p=c2/v**gma
    step3p = p
    step3v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 4-1
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v4=v1
    p4=c2/v4**gma
    p=linspace(p4,p1,NBprocess)
    v=NBprocess*[v1]
    step4p = p
    step4v = v
    
    return step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v


#~~~~~~~~~~~~~~~~~~~
#   Diesel Cycle
#~~~~~~~~~~~~~~~~~~~
def diesel(p_min,p_max,v_max,rc,gma):
    '''This function prints Diesel cycle
    arguments are as follows:
    p_min: minimum pressure
    p_max: Maximum pressure
    v_max: Maximum volume
    rc: Cut-Off ratio
    gma: Adiabatic exponent
    The order of arguments is:
    p_min,p_max,v_max,rc,gma
    
    NOTE: After calling the function you have to type:
    show()
    '''

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 1-2
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p1=p_min
    v1=v_max
    p2=p_max
    v2=v1*(p1/p2)**(1/gma)
    c1=p1*v1**gma
    v=linspace(v1,v2,NBprocess)
    p=c1/v**gma
    step1p = p
    step1v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 2-3
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p3=p2
    p=zeros(NBprocess)
    p=p+p2
    v3=rc*v2
    v=linspace(v2,v3,NBprocess)
    step2p = p
    step2v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 3-4
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v4=v1
    c2=p3*v3**gma
    v=linspace(v3,v4,NBprocess)
    p=c2/v**gma
    step3p = p
    step3v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 4-1
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v4=v1
    v=NBprocess*[v4]
    p4=c2/v4**gma
    p=linspace(p4,p1,NBprocess)
    step4p = p
    step4v = v
    
    return step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v


#~~~~~~~~~~~~~~~~~~~
#    Dual Cycle
#~~~~~~~~~~~~~~~~~~~
def dual(p_min,p_max,v_max,r,rc,gma):
    '''This function prints Dual cycle
    arguments are as follows:
    p_min: minimum pressure
    p_max: Maximum pressure
    v_max: Maximum volume
    r: compression ratio
    rc: Cut-Off ratio
    gma: Adiabatic exponent
    The order of arguments is:
    p_min,p_max,v_max,r,rc,gma
    
    NOTE: After calling the function you have to type:
    show()
    '''
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 1-2
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p1=p_min
    v1=v_max
    c1=p1*v1**gma
    v2=v1/r
    v=linspace(v1,v2,NBprocess)
    p=c1/v**gma
    step1p = p
    step1v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 2-3
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p3=p_max
    p2=c1/v2**gma
    v3=v2
    p=linspace(p2,p3,NBprocess)
    v=zeros(NBprocess)
    v=v+v2
    step2p = p
    step2v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 3-4
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v4=v3*rc
    p4=p3
    p=zeros(NBprocess)
    p=p+p4
    v=linspace(v3,v4,NBprocess)
    step3p = p
    step3v = v
   
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 4-5
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v5=v1
    c2=p4*v4**gma
    p5=c2/v5**gma
    v=linspace(v4,v5,NBprocess)
    p=c2/v**gma
    step4p = p
    step4v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 5-1
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p=linspace(p1,p5,NBprocess)
    v=zeros(NBprocess)
    v=v+v1 
    step5p = p
    step5v = v
   
    return step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v


#~~~~~~~~~~~~~~~~~~~
#    Lenoir Cycle
#~~~~~~~~~~~~~~~~~~~
def lenoir(p_min,p_max,v_max,gma):
    '''This function prints Lenoir cycle
    arguments are as follows:
    p_min: minimum pressure
    p_max: Maximum pressure
    v_max: Maximum volume
    gma: Adiabatic exponent
    The order of arguments is:
    p_min,p_max,v_max,gma
    
    NOTE: After calling the function you have to type:
    show()
    '''
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 3-1
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v1=v_max
    p1=p_min
    p3=p_max
    c=p1*v1**gma
    v3=(c/p3)**(1/gma)
    v=linspace(v1,v3,NBprocess)
    p=c/v**gma
    step1p = p
    step1v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 3-2
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p2=p1
    v2=v3
    p=linspace(p3,p2,NBprocess)
    v=zeros(NBprocess)
    v=v+v2
    step2p = p
    step2v = v


    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 2-1
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v=linspace(v2,v1,NBprocess)
    p=zeros(NBprocess)
    p=p+p2
 
    step3p = p
    step3v = v
    
    return step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v

   

#~~~~~~~~~~~~~~~~~~~
#    Atkinson Cycle
#~~~~~~~~~~~~~~~~~~~
def atkinson(p_min,p_max,v_max,r,gma):
    '''This function prints Atkinson cycle
    arguments are as follows:
    p_min: minimum pressure
    p_max: Maximum pressure
    v_max: Maximum volume (minimum volume just befor compression)
    r: compression ratio
    gma: Adiabatic exponent
    The order of arguments is:
    p_min,p_max,v_max,r,gma
    
    NOTE: After calling the function you have to type:
    show()
    '''

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 1-2
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p1=p_min
    v1=v_max
    c1=p1*v1**gma
    v2=v1/r
    v=linspace(v2,v1,NBprocess)
    p=c1/v**gma
 
    step1p = p
    step1v = v


    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 2-3
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v3=v2
    p3=p_max
    p2=c1/v2**gma
    p=linspace(p2,p3,NBprocess)
    v=zeros(NBprocess)
    v=v+v3

    step1p = p
    step1v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 3-4
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    c2=p3*v3**gma
    p4=p1
    v4=(c2/p4)**(1/gma)
    v=linspace(v3,v4,NBprocess)
    p=c2/v**gma
 
    step2p = p
    step2v = v

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 4-1
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v=linspace(v4,v1,NBprocess)
    p=zeros(NBprocess)+p1 
    step1p = p
    step1v = v
   
    return step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v


#~~~~~~~~~~~~~~~~~~~~~
#    Stirling Cycle
#~~~~~~~~~~~~~~~~~~~~~

def stirling(p_min,p_max,v_max,r,gma):
    '''This function prints Stirling cycle
    arguments are as follows:
    p_min: minimum pressure
    p_max: Maximum pressure
    v_max: Maximum volume 
    r: compression ratio
    gma: Adiabatic exponent
    The order of arguments is:
    p_min,p_max,v_max,r,gma
    
    NOTE: After calling the function you have to type:
    show()
    '''

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 1-2
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p1=p_min
    v1=v_max
    c1=p1*v1
    v2=v1/r
    p2=c1/v2
    v=linspace(v1,v2,NBprocess)
    p=c1/v
    step1p = p
    step1v = v
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 2-3
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p3=p_max
    v3=v2
    v=zeros(NBprocess)+v3
    p=linspace(p2,p3,NBprocess)
    step2p = p
    step2v = v 

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 3-4
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    c2=p3*v3
    v4=v1
    p4=c2/v4
    v=linspace(v3,v4,NBprocess)
    p=c2/v
    step3p = p
    step3v = v
   

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 4-1
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p=linspace(p4,p1,NBprocess)
    v=zeros(NBprocess)+v1
    step4p = p
    step4v = v
   
    return step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v

#~~~~~~~~~~~~~~~~~~~~~
#    Ericsson Cycle
#~~~~~~~~~~~~~~~~~~~~~
def ericsson(p_min,p_max,v_max,r,gma):
    '''This function prints Ericsson cycle
    arguments are as follows:
    p_min: minimum pressure
    p_max: Maximum pressure
    v_max: Maximum volume 
    r: compression ratio
    gma: Adiabatic exponent
    The order of arguments is:
    p_min,p_max,v_max,r,gma
    
    NOTE: After calling the function you have to type:
    show()
    '''

    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 1-2
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p1=p_min
    p2=p1
    v1=v_max
    v2=v1/r
    p=zeros(NBprocess)+p2
    v=linspace(v2,v1,NBprocess)
    step1p = p
    step1v = v
    
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 2-3
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    c1=p2*v2
    p3=p_max
    v3=c1/p3
    v=linspace(v3,v2,NBprocess)
    p=c1/v
    step2p = p
    step2v = v
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 3-4
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    p4=p3
    c2=p1*v1
    v4=c2/p4
    p=zeros(NBprocess)+p4
    v=linspace(v3,v4,NBprocess)
    step3p = p
    step3v = v
 
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    # Process 4-1
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    v=linspace(v4,v1,NBprocess)
    p=c2/v
    step4p = p
    step4v = v
 
    
if model == 'Stirling':
    step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v = stirling(1e5,20e5,0.5,compression,kappa)    
if model == 'Carnot':
    step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v = carnot(1e5,20e5,0.5,compression,kappa)
if model == 'Otto':
    step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v = otto(1e5,20e5,0.5,compression,kappa)
if model == 'Diesel':
    step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v = diesel(1e5,18e5,0.5,compression,kappa)


NkB = step1p[0]*step1v[0]/300    

fig, ax = plt.subplots(1,2,figsize=(9,4.5))
fig.suptitle(model + '-cycle')
ax[0].set(xlabel='volume (m$^3$)')
ax[0].set(ylabel='pressure (bar)')
ax[0].set(xlim = (0,0.6),ylim=(0,22))

if process1 == 'adiabatic': colorp = 'r'
if process1 == 'isotherm': colorp = 'orange'
if process1 == 'isochor': colorp = 'blue'
if process1 == 'isobar': colorp = 'g'
line1 = ax[0].plot(0,0,lw=2,color=colorp,label=process1)[0]

if process2 == 'adiabatic': colorp = 'r'
if process2 == 'isotherm': colorp = 'orange'
if process2 == 'isochor': colorp = 'blue'
if process2 == 'isobar': colorp = 'g'
line2 = ax[0].plot(0,0,lw=2,color=colorp,label=process2)[0]

if process3 == 'adiabatic': colorp = 'r'
if process3 == 'isotherm': colorp = 'orange'
if process3 == 'isochor': colorp = 'blue'
if process3 == 'isobar': colorp = 'g'
line3 = ax[0].plot(0,0,lw=2,color=colorp,label=process3)[0]

if process4 == 'adiabatic': colorp = 'r'
if process4 == 'isotherm': colorp = 'orange'
if process4 == 'isochor': colorp = 'blue'
if process4 == 'isobar': colorp = 'g'
line4 = ax[0].plot(0,0,lw=2,color=colorp,label=process4)[0]

ax[0].legend()

text1 = ax[0].text(0.02,0.1,'',fontsize=8,transform=ax[0].transAxes)
textT = ax[0].text(0.02,0.02,'T',fontsize=12,transform=ax[0].transAxes)

# setup info box
infobox = ''
infobox += 'kappa: ' + "{:.2f}".format(kappa) + '\n' 
infobox += 'compression: ' + "{:.0f}".format(compression)
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax[0].text(0.02,0.98,infobox, fontsize=8,bbox=props,verticalalignment='top',transform=ax[0].transAxes)

ax[1].axis('off')
ax[1].set(ylim=(0,0.6),xlim=(-5,5))
patch = ax[1].add_patch(Rectangle((-1,0),2,0.5,ec='b',fc=([0,0,0])))
piston = ax[1].add_patch(Rectangle((-1,0),2,0.01,ec='black',fc='black'))
piston1 = ax[1].add_patch(Rectangle((-0.05,0),0.1,0.1,ec='black',fc='black'))
   
ax[1].add_patch(Rectangle((-1.2,0),0.2,0.6,ec='lightgrey',fc='grey'))
ax[1].add_patch(Rectangle((1,0),0.2,0.6,ec='lightgrey',fc='grey'))

cmap = cm.get_cmap('plasma')

print('Start Simulation')

def animate(k):
    
   global step1p,step1v,step2p,step2v,step3p,step3v,step4p,step4v 
   global patch,piston,piston1
   print(k) 
   #print(step1v)
   #print(step1p/1e5)
   if int(k/NBprocess) == 0:
     line1.set_xdata(step1v[:k])
     line1.set_ydata(step1p[:k]/1e5)
     text1.set_text(process1)
     Temp = step1v[k]*step1p[k]/NkB
     textT.set_text('T (K): '+"{:.0f}".format(Temp))
   
     colornorm = Temp/Tscale+0.2
     if colornorm >1 : colornorm=1 
     colort = cmap(colornorm)# choose a nonlinear color scale  
     patch.remove()
     piston.remove()
     piston1.remove()
     pr = step1v[k]
     patch = ax[1].add_patch(Rectangle((-1,0),2,pr,ec='black',fc=colort))
     piston = ax[1].add_patch(Rectangle((-1,pr),2,0.01,ec='black',fc='black'))
     piston1 = ax[1].add_patch(Rectangle((-0.05,pr),0.1,0.1,ec='black',fc='black'))

   
   if int(k/NBprocess) == 1:
     line2.set_xdata(step2v[:k-NBprocess])
     line2.set_ydata(step2p[:k-NBprocess]/1e5)
     text1.set_text(process2)
     Temp = step2v[k-NBprocess]*step2p[k-NBprocess]/NkB
     textT.set_text('T (K): '+"{:.0f}".format(Temp))

     colornorm = Temp/Tscale+0.2
     if colornorm >1 : colornorm=1 # choose a nonlinear color scale  
     colort = cmap(colornorm)# choose a nonlinear color scale  
     
     patch.remove()
     piston.remove()
     piston1.remove()
     pr = step2v[k-NBprocess]
     patch = ax[1].add_patch(Rectangle((-1,0),2,pr,ec='black',fc=colort))
     piston = ax[1].add_patch(Rectangle((-1,pr),2,0.01,ec='black',fc='black'))
     piston1 = ax[1].add_patch(Rectangle((-0.05,pr),0.1,0.1,ec='black',fc='black'))

 
   if int(k/NBprocess) == 2:
     line3.set_xdata(step3v[:k-2*NBprocess])
     line3.set_ydata(step3p[:k-2*NBprocess]/1e5)
     text1.set_text(process3)
     Temp = step3v[k-2*NBprocess]*step3p[k-2*NBprocess]/NkB
     textT.set_text('T (K): '+"{:.0f}".format(Temp))

     colornorm = Temp/Tscale+0.2
     if colornorm >1 : colornorm=1 # choose a nonlinear color scale  
     colort = cmap(colornorm)# choose a nonlinear color scale  
    
     patch.remove()
     piston.remove()
     piston1.remove()
     pr = step3v[k-2*NBprocess]
     patch = ax[1].add_patch(Rectangle((-1,0),2,pr,ec='black',fc=colort))
     piston = ax[1].add_patch(Rectangle((-1,pr),2,0.01,ec='black',fc='black'))
     piston1 = ax[1].add_patch(Rectangle((-0.05,pr),0.1,0.1,ec='black',fc='black'))

    
   if int(k/NBprocess) == 3:
     line4.set_xdata(step4v[:k-3*NBprocess])
     line4.set_ydata(step4p[:k-3*NBprocess]/1e5)
     text1.set_text(process4)
     Temp = step4v[k-3*NBprocess]*step4p[k-3*NBprocess]/NkB
     textT.set_text('T (K): '+"{:.0f}".format(Temp))
    
     colornorm = Temp/Tscale+0.2
     if colornorm >1 : colornorm=1 # choose a nonlinear color scale 
     colort = cmap(colornorm)# choose a nonlinear color scale  
     
     patch.remove()
     piston.remove()
     piston1.remove()
     pr = step4v[k-3*NBprocess]
     patch = ax[1].add_patch(Rectangle((-1,0),2,pr,ec='black',fc=colort))
     piston = ax[1].add_patch(Rectangle((-1,pr),2,0.01,ec='black',fc='black'))
     piston1 = ax[1].add_patch(Rectangle((-0.05,pr),0.1,0.1,ec='black',fc='black'))

    
   
anim = animation.FuncAnimation(fig,animate,interval=1,frames=NBsteps)
anim.save(animationfile,fps=25,dpi=300)


#diesel(1e5,20e5,0.5,5,1.4)