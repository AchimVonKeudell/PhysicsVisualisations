#
# Etching of a solid
# Examples for plasma pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

model = 'HonSi'
model = 'ConSi'

a0 = 0.0529*1e-9
el = 1.6e-19
amu = 1.6e-27
lamb = 1.309
m = 0.333
q = 0.667

NBParticles = 500
elev = 10

if model == 'HonSi':
    animationfile = 'piper_HonSi.mp4'
    textcollision = 'H on Si'
    m1 = 1
    m2 = 28
    z1 = 1
    z2 = 14
    energy = 10000
if model == 'ConSi':
    animationfile = 'piper_ConSi.mp4'
    textcollision = 'C on Si'
    m1 = 12
    m2 = 28
    z1 = 6
    z2 = 14
    energy = 10000

ecutoff = 3
edisplacement = 3

def f_func(t):
  return lamb*t**(0.5-m)*(1+(2*lamb*t**(1-m))**q)**(-1/q)
              
def t_func(th,eps):
  return (eps*np.sin(th/2))**2

def ds_func(atf,t):
  f = f_func(t)
  return np.pi*(atf)**2*0.5*t**(-1.5)*f

def eps_func(atf,z1,z2,m1,m2,e):
  return atf*m2/(z1*z2*(el)**2*(m1+m2))*e*(4*np.pi*8.85e-12)

def atf_func(z1,z2):
  return a0*0.88534*(np.sqrt(z1)+np.sqrt(z2))**(-2/3)

def thlab_func(m1,m2,th):
  return np.arccos((m1+m2*np.cos(th))/np.sqrt((m1)**2+2*m1*m2*np.cos(th)+(m2)**2));


def kelec_func(z1,z2,m1,m2):
  return 3.83*z1**(7/6)*z2/((z1**(2/3)+z2**(2/3))**(3/2)*m1**(1/2))

# ------------------------------------------------------------------------
# Setup plot
# ------------------------------------------------------------------------
xylim = 50
zlim = 100
xt = []
yt = []
zt = []
fig = plt.figure(figsize=(8,4.5))
axp = fig.add_subplot(1,2,1,projection='3d')
axp.set(xlabel="x",ylabel="y",zlabel="z")
axp.set(xlim=(-xylim,xylim),ylim=(-xylim,xylim),zlim=(0,zlim))
axp.quiver(0,0,0,0,0,1,pivot = 'tip', color="r", length=20,label='impact',alpha=0.5)
axp.text2D(0.05,0.95,textcollision, fontsize=10,transform=axp.transAxes)

# trajectory line
xp = []
yp = []
zp = []
line1 = axp.plot3D(xt,yt,zt,color='r',label='trajectory',lw=0.5)[0]
line1l = axp.plot3D(xt,yt,zt,color='orange',lw=2)[0]
line1p = axp.plot3D(xt,yt,zt,color='b',marker='o',markersize=1,label='stopped ion',lw=0)[0]
axp.legend(fontsize=6)


NBzscale = 100
zscale = np.linspace(0,zlim,NBzscale)
zbin = np.zeros(100)
ax2 = fig.add_subplot(1,2,2)
line2 = ax2.plot(zscale,zbin,color='g',lw=2)[0]
ax2.set(xlabel="depth z",ylabel="p(z)")
ax2.set(xlim=(0,zlim),ylim=(0,0.1))

# info box
infobox = ''
infobox += 'm projectile: ' + "{:.0f}".format(m1) + ' (amu)\n' 
infobox += 'm target: ' + "{:.0f}".format(m2) + ' (amu)\n' 
infobox += 'E_projectil: ' + "{:.0f}".format(energy/1e3) + ' (keV)\n' 
infobox += 'E_cutoff: ' + "{:.0f}".format(ecutoff) + ' (eV)' 
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax2.text(0.05,0.95,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax2.transAxes)


fig.tight_layout(pad=5)

# marker at the end
#xp = np.array(1)
#yp = np.array(1)
#zp = np.array(1)
#line1p = axp.plot3D(xp,yp,zp,color='r',marker='o',markersize=3,lw=0)[0]



# ------------------------------------------------------------------------
# Initialize
# ------------------------------------------------------------------------
atf = atf_func(z1,z2)
n = 1e23*1e6;
stotal = n**(-2/3)
dth = 0.0001
kel = kelec_func(z1,z2,m1,m2)
  

def direction(vx0,vy0,vz0,t,ph):
      sine = np.sqrt(vy0*vy0+vz0*vz0)
      srat = np.sin(t)/sine;
      cx2 = np.cos(t)*vx0+np.sin(t)*sine*np.cos(ph);
      cy2 = np.cos(t)*vy0-srat*(vy0*vx0*np.cos(ph)-vz0*np.sin(ph))
      cz2 = np.cos(t)*vz0-srat*(vy0*vx0*np.cos(ph)-vz0*np.sin(ph))
      un = 1/np.sqrt(cx2*cx2+cy2*cy2+cz2*cz2)
      vx1 = cx2*un
      vy1 = cy2*un
      vz1 = cz2*un+1e-12
      return [vx1,vy1,vz1] 


def animate(k):
  s = 0 
  # start conditions new trajectory
  en = energy*el
  x = []
  x.append(0)
  y = []
  y.append(0)
  z = []
  z.append(0)
  
  print('Particle: ',k,' of ',NBParticles)
  c = 0
  vx = 0
  vy = 0
  vz = 1

  # Follow Trajectory and store in x,y,z
  # until energy is lost
  while en>ecutoff*el:
    c += 1   

    # ----------------------------------------------------
    # nuclear stopping select scattering angle according to
    # cross section
    # ---------------------------------------------------
    r1 = np.random.rand()*0.97
    th = np.pi/2;
    s =0 
    while s<=stotal:
      th = th-dth;
      eps = eps_func(atf,z1,z2,m1,m2,en)
      t = t_func(th,eps)
      if eps*eps-t>=0:
         s = s+ds_func(atf,t)*np.sqrt(t)*np.sqrt(eps*eps-t)*dth;
      if th == 0: break
    
    thmin = th
    sgesamt = s

    th = thmin
    s = 0
    while s/sgesamt<r1:
      th = th + dth
      eps = eps_func(atf,z1,z2,m1,m2,en);
      t = t_func(th,eps);
      if eps*eps-t>=0: 
        s = s+ds_func(atf,t)*np.sqrt(t)*np.sqrt(eps*eps-t)*dth
   

    if th>np.pi/2: th = 0

    deltae = 4*(m1*m2)/(m1+m2)**2*en*(np.sin(th/2))**2;
    etransfer = deltae/el

    # --------------------------------
    #  electronic stopping
    #  -------------------------------
    deltae += kel*np.sqrt(en/el/1e3)*n**(-1/3)*1e15*el*1e-4
    en -=  deltae

    thl = thlab_func(m1,m2,th)
    phi = 2*np.pi*np.random.rand();

    vxn = vx
    vyn = vy
    vzn = vz

    # new direction
    [vx,vy,vz] = direction(vxn,vyn,vzn,thl,phi)

    x.append(x[c-1]+vx)
    y.append(y[c-1]+vy)
    z.append(z[c-1]+vz)
   
    #if z[c]<=0: break

  zbinid = int(z[c-1]/zlim*NBzscale)
  # print(z[c-1])
  if zbinid<NBzscale-1:
      zbin[zbinid] += 1/NBParticles
  line2.set_ydata(zbin)      
  
  # accumulat complete trajectory
  for i in range(len(x)):
      xt.append(x[i])
      yt.append(y[i])
      zt.append(z[i])               
  for i in range(len(x)-1,-1,-1):
      xt.append(x[i])
      yt.append(y[i])
      zt.append(z[i])               
  line1.set_xdata(xt)
  line1.set_ydata(yt)
  line1.set_3d_properties(zt) 
  
  # accumulate last trajectory
  xl = []
  yl = []
  zl = []
  for i in range(len(x)):
      xl.append(x[i])
      yl.append(y[i])
      zl.append(z[i])               
  line1l.set_xdata(xl)
  line1l.set_ydata(yl)
  line1l.set_3d_properties(zl) 
  
  # set final particle
  xp.append(x[c-1])
  yp.append(y[c-1])
  zp.append(z[c-1])               
  line1p.set_xdata(xp)
  line1p.set_ydata(yp)
  line1p.set_3d_properties(zp) 
  
  
  axp.view_init(elev+elev*np.sin(k/NBParticles*2*np.pi),45+k/NBParticles*360)  
            
  fig.suptitle('Particles: ' + "{:.0f}".format(k))

anim = animation.FuncAnimation(fig,animate,interval=1,frames=NBParticles)
anim.save(animationfile,fps=25,dpi=180)

