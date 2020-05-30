import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import pandas as pd

'''
# First set up the figure, the axis, and the plot element we want to animate
fig, ax = plt.subplots(ncols=1, nrows=1)
ax = plt.axes(projection='3d')
line, = ax.plot([], [], lw=2)
x = np.linspace(0, 2, 1000)
y = np.linspace(0, 2, 1000)
z = np.sin(x*y)
ax.plot(x,y,z)

def animate(i):
    line.set_data(x, y)
    line.set_3d_properties(np.sin(x[i]*y[i]))
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate,
                               frames=50, interval=20, blit=True)


plt.show()
'''


filename = r"C:\Users\Salvo\Desktop\IFISC\data\no_words_per_county.csv"
df = pd.read_csv(filename, index_col=[0])

print(df)

df.iloc[1:,1] = df.iloc[1:,1] -1

pd.to_numeric(df.iloc[:,1])
df = df.sort_values(by=['counting'], ascending=False)
df = df.reset_index()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
a = len(df.iloc[:, 1])
x = np.arange(1, a)
#ax.grid(color='xkcd:ocean blue', linestyle=':')
ax.set_yticks(np.arange(0, 71000, 5000), minor=False)
ax.set_xlabel('County ranking')
ax.set_ylabel('No. of nonzero occurrency words')

ax.scatter(x, df.iloc[:-1, 2], s=1)


print(df)

#ax.set_yscale('log')
ax.set_xscale('log')

fig.tight_layout()

plt.show()


#filename = r"C:\Users\Salvo\Desktop\IFISC\data\no_words_per_county(sorted).csv"
#df.to_csv(filename)


'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
line, = ax.plot([], [], [], "o", markersize=2)
X = np.arange(-5, 5, 0.25)
Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)
R = np.sqrt(X**2 + Y**2)
Z = np.sin(R)
i = int
def update(i, x, y, z):
    line.set_data(X[i], Y[i])
    line.set_3d_properties(Z[i])

ani = animation.FuncAnimation(fig, update(i,X,Y,Z), 100, blit = False,)
plt.show()
'''