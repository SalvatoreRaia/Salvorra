import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure()
fig.set_dpi(100)
ax1 = fig.add_subplot(1, 1, 1)

# Speed of the wave
c = 20
scale = 2

# Initial conditions
x0 = np.linspace(0, 3 * np.pi, 1000)
t0 = 0

# Increment
dt = 0.01


# Onda
def u(x, t):
    return scale * np.sin(x - c * t)


a = []

for i in range(500):
    value = u(x0, t0)
    t0 = t0 + dt
    a.append(value)

k = 0


def animate(i):
    global k
    x = a[k]
    k += 1
    ax1.clear()
    plt.plot(x0, x)
    plt.plot(x0, np.ones((x0.shape[0], 1)) * scale)
    plt.plot(x0, np.ones((x0.shape[0], 1)) * (-scale))
    plt.ylim([-scale - 1, scale + 1])


anim = animation.FuncAnimation(fig,animate,frames=360,interval=20)


plt.show()