import numpy as np
import vpython as vp
import matplotlib.pyplot as plt


def VerletUnderground(x_0, v_0, a, t_0, t_1, dt, verlet = 'velocity', want_v = False): # verlet = 'velocity', 'regular'
    X = np.zeros((int((t_1 - t_0)//dt + 1), x_0.shape[0]))
    V = np.zeros((int((t_1 - t_0)//dt + 1), v_0.shape[0]))
    X[0] = x_0
    V[0] = v_0
    if verlet == 'velocity':
        for i in range(1, X.shape[0]):
            X[i] = X[i - 1] + V[i - 1] * dt + 0.5 * a(X[i - 1]) * dt ** 2
            V[i] = V[i - 1] + (a(X[i - 1]) + a(X[i])) * dt / 2
    elif verlet == 'regular':
        X[1] = X[0] + V[0] * dt + 0.5 * a(X[0], V[0]) * dt ** 2
        for i in range(2, X.shape[0]):
            X[i] = 2 * X[i - 1] - X[i - 2] + a(X[i - 1], V[i - 1]) * dt ** 2
    if want_v and verlet == 'velocity':
        return X, V
    else:
        return X


#G = 6.67408*1e-11
#M = 1990000*1e24
#m = 5.974*1e24
#r_0 = np.array([149.6*1e9, 0, 0])
#v_0 = np.array([0, 30.2865*1e3, 0])
r_0 = np.array([1, 0, 0])
v_0 = np.array([0.5, 0.5, 0])
X, V = VerletUnderground(x_0=r_0, v_0=v_0, a=lambda x: -x/np.linalg.norm(x)**3, t_0=0, t_1=100, dt=0.001, want_v=True)
'''
figa, ax = plt.subplots(nrows=1, ncols=1, figsize=(19.20, 10.80))
ax.plot(X.T[0], X.T[1], label='X', color='xkcd:black')
ax.plot(V.T[0], V.T[1], label='V', color='xkcd:wine red')
ax.legend()
#plt.savefig('./potenziale_g_r-3.png')
plt.show()
'''
vp.scene.caption = 'Keplerian Motion'
# vp.scene.forward = vp.vector(0,-.3,-1)

sun = vp.sphere(pos=vp.vector(0, 0, 0), radius=0.1, color=vp.color.yellow)
sun.p = 0
earth = vp.sphere(pos=vp.vector(X[0][0], X[0][1], X[0][2]), radius=0.05, color=vp.color.green, make_trail=True, interval=10, retain=100)
earth.p = vp.vector(V[0][0], V[0][1], V[0][2])
i=0
while True:
    vp.rate(1000)
    sun.p = vp.vector(0, 0, 0)
    earth.p =vp.vector(V[i%X.shape[0]][0], V[i%X.shape[0]][1], V[i%X.shape[0]][2])
    sun.pos = vp.vector(0, 0, 0)
    earth.pos = vp.vector(X[i%X.shape[0]][0], X[i%X.shape[0]][1], X[i%X.shape[0]][2])
    i += 1