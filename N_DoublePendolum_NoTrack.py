

from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

G = 9.81  # acceleration due to gravity, in m/s^2
L1 = 1.0  # length of pendulum 1 in m
L2 = 1.0  # length of pendulum 2 in m
M1 = 1.0  # mass of pendulum 1 in kg
M2 = 1.0  # mass of pendulum 2 in kg


def derivs(state, t):
  
    dydx = np.zeros_like(state)
    dydx[0] = state[1]

    delta = state[2] - state[0]
    den1 = (M1+M2) * L1 - M2 * L1 * cos(delta) * cos(delta)
    dydx[1] = ((M2 * L1 * state[1] * state[1] * sin(delta) * cos(delta)
                + M2 * G * sin(state[2]) * cos(delta)
                + M2 * L2 * state[3] * state[3] * sin(delta)
                - (M1+M2) * G * sin(state[0]))
               / den1)

    dydx[2] = state[3]

    den2 = (L2/L1) * den1
    dydx[3] = ((- M2 * L2 * state[3] * state[3] * sin(delta) * cos(delta)
                + (M1+M2) * G * sin(state[0]) * cos(delta)
                - (M1+M2) * L1 * state[1] * state[1] * sin(delta)
                - (M1+M2) * G * sin(state[2]))
               / den2)

    return dydx

# create a time array sampled at 0.05 second steps
dt = 0.05
t = np.arange(0, 30, dt)

# th1 and th2 are the initial angles (degrees)
# w10 and w20 are the initial angular velocities (degrees per second)
th1 = 160.0001
w1 = 0.0
th2 = -10.0001
w2 = 0

N = 40 # N pendulums

x = [  ]
y = [  ]

for i in range (N):

    # initial state
    displacement = 0.001
    state = np.radians([th1, w1, th2 + displacement*i, w2])

    # integrate your ODE using scipy.integrate.
    sol = integrate.odeint(derivs, state, t)

    x1 = L1*sin(sol[:, 0])
    y1 = -L1*cos(sol[:, 0])

    x2 = L2*sin(sol[:, 2]) + x1
    y2 = -L2*cos(sol[:, 2]) + y1

    x.append([ x1 , x2 ])
    y.append([ y1 , y2 ])


fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax.set_aspect('equal')
ax.set_title('Wait For it')
ax.grid()

#line, = ax.plot([], [], 'o-', lw=2)

#lines = sum([ax.plot([], [], 'o-')], [], lw=2)

varacaso = []
for i in range(N):
    varacaso.append(ax.plot([], [], 'o-'))

lines = varacaso   

for i in range(N):
    varacaso.append(ax.plot([], [], '--'))
pts = varacaso  

time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

track_segments = np.zeros((N, 0, 2))
tracks = collections.LineCollection(N)
tracks.set_array(np.linspace(0, 1, N))
ax.add_collection(tracks)


NN = range(N)


def init():

    for line, ptss in zip (lines, pts) :
        #print(line[0])
        line[0].set_data([], [])
        #ptss[0].set_data([], [])
    
    tracks.set_segments(np.zeros((N, 0, 2)))

    time_text.set_text('')
    return lines, time_text, tracks #, pts


path = 10

def animate(i):
    if i > path :
        for n, line, ptss in zip (NN, lines, pts):
            thisx = [0, x[n][0][i], x[n][1][i-path:i]]
            thisy = [0, y[n][0][i], y[n][1][i-path:i]]

            #ptss[0].set_data(thisx, thisy)
            line[0].set_data(thisx[-1], thisy[-1])

    else :

        for n, line, ptss in zip (NN, lines, pts):
            thisx = [0, x[n][0][i], x[n][1][i]]
            thisy = [0, y[n][0][i], y[n][1][i]]

            line[0].set_data(thisx, thisy)
            #ptss[0].set_data(thisx, thisy)

    sl = slice(max(0, i - track_length), i)
    tracks.set_segments(positions[:, sl, -1])

    time_text.set_text(time_template % (i*dt))

    return lines, time_text, tracks #pts


ani = animation.FuncAnimation(fig, animate, range(1, len(sol)),
                              interval=dt*0.1, blit=False)#, init_func=init)
plt.show()

"""

"""