#
# Double pendulum
# Trajectory from Hamiltonian solution to the problem
# see for example:
# https://scienceworld.wolfram.com/physics/DoublePendulum.html    
# Examples for pyhsics lectures
# Achim von Keudell
# Ruhr University Bochum, 2022
#

import matplotlib.pyplot as plt 
import numpy as np
import matplotlib.animation as animation

model = 'equal pendulum'  
model = 'unequal pendulum m2 large'  
model = 'unequal pendulum m1 large'  
       
l = 1
g = 9.8  
dtplot = 0.01
tend = 200
steps = int(tend/dtplot)

# time step has to be set rather small 
# so that numerical error do not accumulate 
# too much, otherwise the phase space 
# will not be covered fully. This is adjusted by
# sublooping = timesteps in between plotting
sublooping = 2000
dt = dtplot/sublooping 
variablemasses = True

if model == 'equal pendulum':
   m1 = 1
   m2 = 1
   l1 = 1
   l2 = 1        
   animationfile = 'doublependulum_equal.mp4'  
if model == 'unequal pendulum m2 large':
   m1 = 1
   m2 = 2
   l1 = 1.5
   l2 = 0.5        
   animationfile = 'doublependulum_m2_large.mp4'      
if model == 'unequal pendulum m1 large':  
   m1 = 2
   m2 = 1
   l1 = 0.5
   l2 = 1.5        
   animationfile = 'doublependulum_m1_large.mp4'  

p1 = 0
p2 = 0
theta1 = np.pi
theta2 = np.pi-0.01
    
def pendulumposition(theta1,p1,theta2,p2):
        x1 =  l * np.sin(theta1)        
        y1 = -l * np.cos(theta1)
          
        x2 = x1 + l * np.sin(theta2)
        y2 = y1 - l * np.cos(theta2)
         
        #print(self.theta1, self.theta2)
        return np.array([[0.0, 0.0], [x1, y1], [x2, y2]])

def pendulumpositionvariable(theta1,p1,theta2,p2):
        x1 =  l1 * np.sin(theta1)        
        y1 = -l1 * np.cos(theta1)
          
        x2 = x1 + l2 * np.sin(theta2)
        y2 = y1 - l2 * np.cos(theta2)
         
        #print(self.theta1, self.theta2)
        return np.array([[0.0, 0.0], [x1, y1], [x2, y2]])
        
     
def pendulumstepvariable():
         
        global p1,p2,theta1,theta2 
  
        expr1 = np.cos(theta1 - theta2)
        expr2 = np.sin(theta1 - theta2)
        expr22 = (m1+m2*expr2**2)
        
        A = p1 * p2 * expr2 / (l1*l2*expr22)
        B = ((l2**2*m2*p1**2 + l1**2*(m1+m2)*p2**2 - l1*l2*m2*p1 * p2 * expr1)* np.sin(2 * (theta1 - theta2)) / 
               (2*l1**2*l2**2*expr22**2))
         
        theta1 += dt * (p1*l2 - l1*p2 * expr1) / (l1**2*l2*expr22)
        theta2 += dt * (l1*(m1+m2)*p2 - l2*m2*p1 * expr1) / (l1*l2**2*m2*expr22)
        p1 += dt * (-(m1+m2) * g * l1 * np.sin(theta1) - A+B)
        p2 += dt * (-m2*g * l2 * np.sin(theta2) + A - B)
         
        new_position = pendulumpositionvariable(theta1,p1,theta2,p2)
        return new_position


def pendulumstep():
         
        global p1,p2,theta1,theta2 

        expr1 = np.cos(theta1 - theta2)
        expr2 = np.sin(theta1 - theta2)
        expr3 = (1 + expr2**2)
        expr4 = p1 * p2 * expr2 / expr3
        expr5 = (p1**2 + 2 * p2**2 - p1 * p2 * expr1) \
        * np.sin(2 * (theta1 - theta2)) / 2 / expr3**2
        expr6 = expr4 - expr5
         
        theta1 += dt * (p1 - p2 * expr1) / expr3
        theta2 += dt * (2 * p2 - p1 * expr1) / expr3
        p1 += dt * (-2 * g * l * np.sin(theta1) - expr6)
        p2 += dt * (    -g * l * np.sin(theta2) + expr6)
         
        new_position = pendulumposition(theta1,p1,theta2,p2)
        return new_position
 
# setup plot

fig, ax = plt.subplots(figsize=(4.5,4.5))
ax.set(xlabel='x (m)',ylabel='y (m)')
ax.set_ylim(-2.5, 2.5)
ax.set_xlim(-2.5, 2.5)
trace = ax.plot(0,0, linestyle='dashed',lw=0.5,color='orange')[0]
pendulum = ax.plot(0,0)[0]
pendulummarker1 = ax.plot(0,0, marker='o',markersize=m1*5,lw=0)[0]
pendulummarker2 = ax.plot(0,0, marker='o',markersize=m2*5,lw=0)[0]

tracex = []
tracey = []

# info box
infobox = ''
infobox += 'm1: '+"{:.1f}".format(m1)+' (kg) \n'
infobox += 'l1: '+"{:.1f}".format(l1)+' (m) \n'
infobox += 'm2: '+"{:.1f}".format(m2)+' (kg) \n'
infobox += 'l2: '+"{:.1f}".format(l2)+' (m)'
#infobox += 'v_phase: '+ "{:.0f}".format(vphase)+' (m/s)'
props = dict(boxstyle='round', facecolor='lightblue', alpha=0.5) 
ax.text(0.02,0.98,infobox, fontsize=6,bbox=props,verticalalignment='top',transform=ax.transAxes)
  
t = 0 
     
def animate(k):
    
    global p1,p2,theta1,theta2 
    global t,steps
    
    print('time: '+ "{:.2f}".format(t)+' s of '+"{:.0f}".format(steps*dtplot)+' s ')
    
    t += dt*sublooping
    
    for i in range(sublooping):
       if variablemasses: 
         new_position = pendulumstepvariable()
       else:
         new_position = pendulumstep()
  
    tracex.append(new_position[2,0])
    tracey.append(new_position[2,1])
    trace.set_xdata(tracex)
    trace.set_ydata(tracey)
    pendulum.set_xdata(new_position[:, 0])
    pendulum.set_ydata(new_position[:, 1])
    
    pendulummarker1.set_xdata(new_position[1, 0])
    pendulummarker1.set_ydata(new_position[1, 1])
    
    pendulummarker2.set_xdata(new_position[2, 0])
    pendulummarker2.set_ydata(new_position[2, 1])
    
    fig.suptitle('Time: ' + "{:.2f}".format(t)+' s')
        
        
anim = animation.FuncAnimation(fig,animate,interval=0,frames=steps)
anim.save(animationfile,fps=25,dpi=300) 
