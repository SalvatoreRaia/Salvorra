import numpy as np
import vpython as vp
import matplotlib.pyplot as plt


def HeunkHeunk(x_0, f, t_0, t_1, dt, f_tuple=None, system=False): # {HEUN METHOD}
     if system: # f = f(t, Y, f_tuple) dove Y è la matrice delle derivate
         X = np.zeros((int((t_1 - t_0) // dt + 1), x_0.shape[0], x_0.shape[1]))
         X[0] = x_0
         for i in range(1, X.shape[0]):
            for j in range(X.shape[1] - 1):
                x_tilde = X[i - 1][j] + dt * X[i - 1][j + 1]
                X[i] = X[i - 1][j] + 0.5 * dt * (X[i - 1][j] + x_tilde)
            x_tilde = X[i - 1] + dt * f(t_0 + dt * (i - 1), X[i - 1], f_tuple)
            X[i] = X[i - 1][j] + 0.5 * dt * (f(t_0 + dt * (i - 1), X[i - 1], f_tuple) + f(t_0 + dt * i, x_tilde, f_tuple))
     else: # f = f(t, y, f_tuple)
        X = np.zeros((int((t_1 - t_0)//dt + 1), x_0.shape[0]))
        X[0] = x_0
        for i in range(1, X.shape[0]):
            x_tilde = X[i-1] + dt*f(t_0+dt*(i-1), X[i-1], f_tuple)
            X[i] = X[i-1] + 0.5*dt*(f(t_0+dt*(i-1), X[i-1], f_tuple) + f(t_0+dt*i, x_tilde, f_tuple))
     return X


def LarmorTau(t, L, f_tuple=(1, np.array([0, 0, 1]))): # f_tuple = (gamma, B) => L' = gamma*L x B
    return np.cross(f_tuple[0]*L, f_tuple[1])


gamma = -1
B = np.array([0, 0, 0.1])
dt = 0.01
L = HeunkHeunk(x_0=np.array([1, 1, 2]), f=LarmorTau, t_0=0, t_1=400, dt=dt, f_tuple=(gamma, B), system=False)


figa, ax = plt.subplots(nrows=1, ncols=1)
ax.plot(L.T[0], label='$L_x$', color='xkcd:blurple')
ax.plot(L.T[1], label='$L_y$', color='xkcd:moss green')
ax.plot(L.T[2], label='$L_z$', color='xkcd:fire engine red')
plt.show()

vp.scene.forward = vp.vector(0, -1, 0)
el_pos = np.array([0, 1., 0])
el_v = np.cross(el_pos, L[0])/np.dot(L[0], L[0])
momentum = vp.arrow(pos=vp.vector(0, 0, 0), axis=vp.vector(L[0][0], L[0][1], L[0][2]), color=vp.color.blue, shaftwidth=0.05)
B_vect = vp.arrow(pos=vp.vector(0, 0, 0), axis=vp.vector(B[0], B[1], B[2]), color=vp.color.red, shaftwidth=0.01)
tau = LarmorTau(t=0, L=L[0], f_tuple=(gamma, B))
momentum.p = vp.vector(tau[0], tau[1], tau[2])
electron = vp.sphere(pos=vp.vector(el_pos[0], el_pos[1], el_pos[2]), radius=0.05, color=vp.color.green, make_trail=True, interval=10, retain=100)
electron.p = vp.vector(el_v[0], el_v[1], el_v[2])
i=0
while True:
    vp.rate(100)
    B_vect = vp.arrow(pos=vp.vector(0, 0, 0), axis=vp.vector(14*B[0], 14*B[1], 14*B[2]), color=vp.color.red, shaftwidth=0.05)
    el_pos += el_v*dt
    el_v = np.cross(el_pos, L[i]) / np.dot(L[i], L[i])
    tau = LarmorTau(t=0, L=L[i], f_tuple=(gamma, B))
    momentum.p = vp.vector(tau[0], tau[1], tau[2])
    momentum.axis = vp.vector(L[i][0], L[i][1], L[i][2])
    electron.pos = vp.vector(el_pos[0], el_pos[1], el_pos[2])
    electron.p = vp.vector(el_v[0], el_v[1], el_v[2])
    i += 1
