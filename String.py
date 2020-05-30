import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
import matplotlib.animation as animation


fig = plt.figure()
fig.set_dpi(100)
ax1 = fig.add_subplot(1, 1, 1)

# Wave speed
c = 1

# x axis
x0 = np.linspace(-pi, pi, 100)

# Initial time
t0 = 0

# Time increment
dt = 0.1


# Wave equation solution
def u(x, t):
    return 0.5 * (np.sin(x + c * t) + np.sin(x - c * t))


a = []

for i in range(100):
    value = u(x0, t0)
    t0 = t0 + dt
    a.append(value)

k = 0

print(a)
def animate(i):
    global k
    x = a[k]
    k += 1
    ax1.clear()
    plt.plot(x0, x, color='cyan')
    plt.grid(True)
    plt.ylim([-2, 2])
    plt.xlim([-pi, pi])


anim = animation.FuncAnimation(fig, animate, frames=120, interval=20, repeat=True)
plt.show()