import scipy as sci
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import numpy as np
from matplotlib.colors import cnames

"""
ANIMATION FUNCTION
   
   Update the data held by the scatter plot and therefore animates it.
   Args:
       iteration (int): Current iteration of the animation
       data (list): List of the data positions at each iteration.
       scatters (list): List of all the scatters (One per element)
   Returns:
       list: List of scatters (One per element) with new coordinates
   """


#Define universal gravitation constant
G=6.67408e-11 #N-m2/kg2
#Reference quantities
m_nd=1.989e+30 #kg #mass of the sun
r_nd=5.326e+12 #m #distance between stars in Alpha Centauri
v_nd=30000 #m/s #relative velocity of earth around the sun
t_nd=79.91*365*24*3600*0.51 #s #orbital period of Alpha Centauri
#Net constants
K1=G*t_nd*m_nd/(r_nd**2*v_nd)
K2=v_nd*t_nd/r_nd

#Define masses
m1=1.1 #Alpha Centauri A
m2=0.907 #Alpha Centauri B
m3=0.9 #Third Star

sizes = [10*m1, 10*m2, 10*m3]

#Define initial position vectors
r1=[-0.5,0,0] #m
r2=[0.5,0,0] #m
r3=[0,1,0] #m

#Convert pos vectors to arrays
r1=sci.array(r1,dtype="float64")
r2=sci.array(r2,dtype="float64")
r3=sci.array(r3,dtype="float64")

#Find Centre of Mass
#r_com=(m1*r1+m2*r2)/(m1+m2)
r_com=(m1*r1+m2*r2+m3*r3)/(m1+m2+m3)

#Define initial velocities
v1=[0.01,0.01,0] #m/s
v2=[-0.05,0,-0.1] #m/s
v3=[0,-0.01,0]

#Convert velocity vectors to arrays
v1=sci.array(v1,dtype="float64")
v2=sci.array(v2,dtype="float64")
v3=sci.array(v3,dtype="float64")

#Find velocity of COM
#v_com=(m1*v1+m2*v2)/(m1+m2)
v_com=(m1*v1+m2*v2+m3*v3)/(m1+m2+m3)


#A function defining the equations of motion 
def TwoBodyEquations(w,t,G,m1,m2):
    r1=w[:3]
    r2=w[3:6]
    v1=w[6:9]
    v2=w[9:12]
    r=sci.linalg.norm(r2-r1) #Calculate magnitude or norm of vector
    dv1bydt=K1*m2*(r2-r1)/r**3
    dv2bydt=K1*m1*(r1-r2)/r**3
    dr1bydt=K2*v1
    dr2bydt=K2*v2
    r_derivs=sci.concatenate((dr1bydt,dr2bydt))
    derivs=sci.concatenate((r_derivs,dv1bydt,dv2bydt))
    return derivs

def ThreeBodyEquations(w,t,G,m1,m2,m3):
    r1=w[:3]
    r2=w[3:6]
    r3=w[6:9]
    v1=w[9:12]
    v2=w[12:15]
    v3=w[15:18]
    r12=sci.linalg.norm(r2-r1)
    r13=sci.linalg.norm(r3-r1)
    r23=sci.linalg.norm(r3-r2)
    
    dv1bydt=K1*m2*(r2-r1)/r12**3+K1*m3*(r3-r1)/r13**3
    dv2bydt=K1*m1*(r1-r2)/r12**3+K1*m3*(r3-r2)/r23**3
    dv3bydt=K1*m1*(r1-r3)/r13**3+K1*m2*(r2-r3)/r23**3
    dr1bydt=K2*v1
    dr2bydt=K2*v2
    dr3bydt=K2*v3
    r12_derivs=sci.concatenate((dr1bydt,dr2bydt))
    r_derivs=sci.concatenate((r12_derivs,dr3bydt))
    v12_derivs=sci.concatenate((dv1bydt,dv2bydt))
    v_derivs=sci.concatenate((v12_derivs,dv3bydt))
    derivs=sci.concatenate((r_derivs,v_derivs))
    return derivs

"""
#Package initial parameters
init_params=sci.array([r1,r2,v1,v2]) #create array of initial params
init_params=init_params.flatten() #flatten array to make it 1D
time_span=sci.linspace(0,8,500) #8 orbital periods and 500 points
#Run the ODE solver
import scipy.integrate
two_body_sol=sci.integrate.odeint(TwoBodyEquations,init_params,time_span,args=(G,m1,m2))

r1_sol=two_body_sol[:,:3]
r2_sol=two_body_sol[:,3:6]
"""

#Package initial parameters
init_params=sci.array([r1,r2,r3,v1,v2,v3]) #Initial parameters
init_params=init_params.flatten() #Flatten to make 1D array
time_span=sci.linspace(0,35,3500) #20 orbital periods and 500 points
#Run the ODE solver

print('suca')

import scipy.integrate
three_body_sol=sci.integrate.odeint(ThreeBodyEquations,init_params,time_span,args=(G,m1,m2,m3))

r1_sol=three_body_sol[:,:3]
r2_sol=three_body_sol[:,3:6]
r3_sol=three_body_sol[:,6:9]

x1 = r1_sol[:,0]
x2 = r2_sol[:,0]
x3 = r3_sol[:,0]
y1 = r1_sol[:,1]
y2 = r2_sol[:,1]
y3 = r3_sol[:,1] 
z1 = r1_sol[:,2]
z2 = r2_sol[:,2]
z3 = r3_sol[:,2]

data = np.array([[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]])

#Create figure
fig=plt.figure()
#Create 3D axes
ax=fig.add_subplot(111,projection="3d")

"""
#Plot the orbits
ax.plot(r1_sol[:,0],r1_sol[:,1],r1_sol[:,2],color="darkblue")
ax.plot(r2_sol[:,0],r2_sol[:,1],r2_sol[:,2],color="tab:red")
ax.plot(r3_sol[:,0],r3_sol[:,1],r3_sol[:,2],color="green")

#Plot the final positions of the stars
ax.scatter(r1_sol[0,0],r1_sol[0,1],r1_sol[0,2],color="darkblue",marker="o",s=100,label="Alpha Centauri A")
ax.scatter(r2_sol[0,0],r2_sol[0,1],r2_sol[0,2],color="tab:red",marker="o",s=100,label="Alpha Centauri B")
ax.scatter(r3_sol[0,0],r3_sol[0,1],r3_sol[0,2],color="green",marker="o",s=100,label="Alpha Centauri C")


"""

#Add a few more bells and whistles
ax.set_xlabel("x",fontsize=14)
ax.set_ylabel("y",fontsize=14)
ax.set_zlabel("z",fontsize=14)
ax.set_title("Visualization of orbits of stars in a 3-body system\n",fontsize=14)
#ax.legend(loc="upper left",fontsize=14)


xmin=x3.min(); xmax=x3.max()
ymin=y3.min(); ymax=y3.max()        
zmin=z3.min(); zmax=z3.max()
ax.set_xlim(xmin-0.1*(xmax-xmin),xmax+0.1*(xmax-xmin))
ax.set_ylim(ymin-0.1*(ymax-ymin),ymax+0.1*(ymax-ymin))
ax.set_zlim(zmin-0.1*(zmax-zmin),zmax+0.1*(zmax-zmin))
#ax.axis('off')
#fig.set_facecolor('black')
#ax.set_facecolor('black')

N_trajectories = 3
colors = plt.cm.jet(np.linspace(0, 1, N_trajectories))

# nOTE: Can't pass empty arrays into 3d version of plot()
#lines = [ax.plot(dat[0, 0:2], dat[1, 0:2], dat[2, 0:2], ls ='-')[0] for dat in data]

#Alternative delcaration of lines objects to be cheked

lines = sum([ax.plot([], [], [], '-', c=c)
             for c in colors], [])


#pts = sum([ax.plot([], [], [], 'o', c=c)
#           for c in colors], [])

varacaso = []
for c, s in zip (colors, sizes):
    varacaso.append(ax.plot([], [], [], 'o', c=c, markersize=s))

pts = sum(varacaso, [])     
    


#time_text = sum([ax.text(0.02, 0.95, 0.8, '')], [])

#ax.view_init(30, 0)

def animate(i, data, lines, pts):

    for line, pt, data in zip(lines, pts, data):
        
        # nOTE: there is no .set_data() for 3 dim data...
        if i > 150 :
            line.set_data(data[0:2, i-150:i])
            line.set_3d_properties(data[2, i-150:i])
        else :
            line.set_data(data[0:2, :i])
            line.set_3d_properties(data[2, :i])
        #line.set_marker("-o")

        pt.set_data(data[0:2, i-2])
        pt.set_3d_properties(data[2, i-2])

        ax.text(0.02, 0.95, 0.8, 'time = %.1fy' % 0.51*i)
        # camera moves
        #ax.view_init(30, 0.3 *i)
        #ax.autoscale_view()
        #fig.canvas.draw()
    #time_text.set_text('time = %.1fy' % 0.51*i)
           
    return  lines + pts#, time_text

ani = animation.FuncAnimation(fig, animate, frames= 1000000,
                               blit=True, fargs = (data, lines, pts), repeat=False, interval = 0.1)#, init_func=init)

#Save as mp4. This requires mplayer or ffmpeg to be installed
#ani.save('lorentz_attractor.mp4', fps=60, extra_args=['-vcodec', 'libx264'])
#ani.save('C:\\Users\\Salvo\\Desktop\\animation2.gif', writer='imagemagick', fps=60)
plt.show()
#Writer = animation.writers['ffmpeg']
#writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800, extra_args=['-vcodec', 'libx264'])
#ani.save('.\Pulla.mp4', writer=writer)

"""

    #perform animation step
    global pendulum, dt
    pendulum.step(dt)
    
    #pendulum.pos ritorna (x,y) con x e y array 1D

    line.set_data(*pendulum.position())
    time_text.set_text('time = %.1f' % pendulum.time_elapsed)
"""
    #l1.set_data(r1_sol[:i,0],r1_sol[:i,1],r1_sol[:i,2])
    #l2.set_data(r2_sol[:i,0],r2_sol[:i,1],r2_sol[:i,2])
    #l3.set_data(r3_sol[:i,0],r3_sol[:i,1],r3_sol[:i,2])
"""
    l1.set_data(x1[:i], y1[:i], z1[:i])
    l2.set_data(x2[:i], y2[:i], z2[:i])
    l3.set_data(x3[:i], y3[:i], z3[:i])
"""

    #time_text.set_text(time_template % (num*dt))
    #return earthDot,venusDot,l1,l2,time_text

       #, time_text

# choose the interval based on dt and the time to animate one step
#from time import time
#t0 = time()
#animate(0)
#t1 = time()
#interval = 1000 * dt - (t1 - t0)



#interval=interval,


"""
t = np.arange(0, 100, dt)
from scipy.integrate import odeint
y = odeint(pend, y0tv, t,atol=1e-16)
y1 = y[:, 1]
x1 = y[:, 0]
y2 = y[:, 3]
x2 = y[:, 2]
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-100, 100),
ylim=(-100, 100))
ax.set_aspect('equal')
#ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.005, 0.9, '', transform=ax.transAxes)
#redDot, = plt.plot([x1[0]], [y1[0]], 'yo')
earthDot, = plt.plot([0], [0], 'bo')
venusDot, = plt.plot([0], [0], 'yo')
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text
l1, = plt.plot([], [], 'r-')
l2, = plt.plot([], [], 'y-')
#l2, =ax.plot(x2, y2, color="yellow")
def animate(num,x1,x2,y1,y2, l1,l2):
    #thisx = [0, x1[i]]
    #thisy = [0, y1[i]]

    #line.set_data(thisx, thisy)
    #thisx = [0, x1[i]]
    #thisy = [0, y1[i]]
    venusDot.set_data(x1[num], y1[num])
    earthDot.set_data(x2[num], y2[num])
    l1.set_data(x1[:num],y1[:num])
    l2.set_data(x2[:num],y2[:num])
    time_text.set_text(time_template % (num*dt))
    return earthDot,venusDot,l1,l2,time_text





    for line, data in zip(lines, dataLines):
        # nOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])

        # Lines to plot in 3D
t = np.linspace(-2*np.pi,2*np.pi,50)
x1, y1, z1 = np.cos(t), np.sin(t), t/t.max()
x2, y2, z2 = t/t.max(), np.cos(t), np.sin(t)
data = np.array([[x1,y1,z1],[x2,y2,z2]])

# NOTE: Can't pass empty arrays into 3d version of plot()
lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]


line_ani = animation.FuncAnimation(fig, update_lines, 50, fargs=(data, lines),
                                   interval=100, blit=True, repeat=True)
"""

#https://pythonmatplotlibtips.blogspot.com/2017/12/draw-3d-line-animation-using-python-matplotlib-funcanimation.html