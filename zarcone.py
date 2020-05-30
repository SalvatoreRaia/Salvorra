import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


x = np.arange(0, 2 * np.pi, 0.01)
fig, ax = plt.subplots()
line = ax.plot(np.sin(x))
arrow = ax.arrow(np.e,0, 0,np.sin(np.e))
ax.set_xlim([0.0, 2 * np.pi])
i=0
x1 = 0
#redDot = plt.plot([x1], [np.sin(0.01)], 'ro')



def update(i, line):
    dt = 0.1
    x = np.arange(0, 2 * np.pi, 0.01)
    A = np.sin(6*x)
    line.set_data(x, A*np.cos(i*dt))
    return line


ani = animation.FuncAnimation(fig, update, 120, fargs=(line), interval=20, blit=False, repeat=True)
plt.show()