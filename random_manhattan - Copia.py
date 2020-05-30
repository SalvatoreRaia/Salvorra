import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import matplotlib.pylab as pl


def inerzia(pindex, p=None, n_dir=4):
    # p probabilità che vada in direzione di inerzia, pindex direzione di inerzia ([x, y, -x, -y] indici antiorari)
    P = np.zeros(n_dir)
    if p is None:
        P = None
    else:
        # assegno equiprobabilità a tutte le altre direzioni
        P[pindex] = p
        for i in range(n_dir):
            if i != pindex:
                P[i] = (1 - p) / (n_dir - 1)
        if np.sum(P) != 1:
            P[pindex] += 1 - np.sum(P)
    return P


def democrazia(array):
    p = np.bincount(array)/array.shape[0]




def stepper(n_dim, step):
    X = np.zeros((n_dim*2, n_dim))
    for i in range(n_dim):
        X[i][i] = step
        X[i+n_dim][i] = -step
    return X

def absolute_stepper(n_dim, step):
    return step*np.identity(n_dim)


def MaxwellBoltzmann(v, m, T): # v_picco = sqrt(2kT/m), v_rms = sqrt(3kT/m)
    kB = 1.38e-23
    return (m/(2*np.pi*kB*T))**1.5 * 4*np.pi * v**2 * np.exp(-m*v**2/(2*kB*T))


def DistributoreDiMaxwell(v_max, res, num_v, T=293.15, m=1.0): #m in uma, T in kelvin
    V = np.arange(0, v_max + res, res)
    MB = np.zeros(V.shape[0])
    VP = np.zeros(num_v)
    for i in range(V.shape[0]):
        MB[i] = MaxwellBoltzmann(v=V[i], T=T, m=m * 1.6605e-27)
    for i in range(1, VP.shape[0] + 1):
        VP[i - 1] = np.trapz(MB[(i - 1) * v_max:i * v_max], V[(i - 1) * v_max:i * v_max])
    scarto = (1 - np.sum(VP)) / VP.shape[0]
    VP += scarto
    print("Ho fatto l'integrale!")
    return VP


class RandomManhattan(object):

    def __init__(self, n_dim=2, spatial=True, p_inerzia=None, n_steps=1000, stepper=stepper): #stepper(n_dim, step)
        self.n_dim = n_dim
        if spatial:
            self.n_dir = n_dim*2
        else:
            self.n_dir = n_dim
        self.p_inerzia = p_inerzia
        self.n_steps = n_steps
        self.stepper = stepper

    def LukeRandomwalker(self, step=lambda: 1, pidx=None, random_seed=time.time()):
        np.random.seed(int(random_seed))
        X = np.zeros((self.n_steps, self.n_dim))
        S = self.stepper(self.n_dim, step())
        if pidx == None:
            pidx = np.random.choice(np.arange(self.n_dir))
        for i in range(1, self.n_steps):
            previous = pidx
            P = inerzia(pindex=pidx, p=self.p_inerzia, n_dir=self.n_dir)
            pidx = np.random.choice(np.arange(self.n_dir), p = P)
            if pidx != previous:
                S = stepper(self.n_dim, step())
            new_step = S[pidx]
            X[i] = X[i - 1] + new_step
        return X, pidx

    def EnderMan(self, N, step=lambda: 1, pidx=None):
        E = np.zeros((N, self.n_dim))
        random_seed = time.time()
        for i in range(N):
            X, _ = self.LukeRandomwalker(step=step, random_seed=random_seed, pidx=pidx)
            E[i] = X[-1]
            random_seed += 7
            print('Ripetizioni effettuate: %d / %d' % (i+1, N))
        return E

    def Elezioni(self, n_citt, odds, step=lambda: 1, random_seed = time.time()):
        np.random.seed(int(random_seed))
        E = np.zeros((n_citt, self.n_dim))
        PIDX = np.zeros(n_citt)
        for i in range(self.n_steps):
            for n in range(n_citt):
                previous = PIDX[n]
                P = inerzia(pindex=PIDX[n], p=self.p_inerzia, n_dir=self.n_dir)
                PIDX[n] = np.random.choice(np.arange(self.n_dir), p = P)
                if PIDX[n] != previous:
                    S = stepper(self.n_dim, step())
                new_step = S[PIDX[n]]
                E[n][i] = E[n][i - 1] + new_step
                print('Faccio propaganda: ', i*n_citt + n + 1, '/', n_steps*n_citt)
            odds = PIDX
            #ricalcolo
        return E, PIDX


#VP = DistributoreDiMaxwell(10000, 0.01, 100)
n_citt = 1000
n_steps = 20
n_dim = 42
myranda = RandomManhattan(n_dim=n_dim, spatial=False, p_inerzia=0.25, n_steps=n_steps, stepper=absolute_stepper)
Popolo = np.zeros((n_citt, n_steps, n_dim))
LastVote = np.zeros(n_citt)
for n in range(Popolo.shape[0]):
    Popolo[n], LastVote[n] = myranda.LukeRandomwalker(step=lambda: 1, random_seed=time.time()+n, pidx=39)
    print('Faccio propaganda: ', n+1, '/', n_citt)
colores = pl.cm.gist_rainbow(np.linspace(0, 1, n_citt))

#X = myranda.LukeRandomwalker(step=lambda: 1, pidx=0)
#Y = myranda.LukeRandomwalker(step=lambda: 1, random_seed=3450954, pidx=1)
#Z = myranda.LukeRandomwalker(step=lambda: 1, random_seed=2057, pidx=2)
'''
figa, ax = plt.subplots(ncols=1, nrows=1)
ax.plot(X.T[0], X.T[1], zorder=1)
ax.scatter(X[-1][0], X[-1][1], color='xkcd:fire engine red', s=np.max(X.flatten())/1000, zorder=100)
ax.plot(Y.T[0], Y.T[1], zorder=1)
ax.scatter(Y[-1][0], Y[-1][1], color='xkcd:fire engine red', s=np.max(X.flatten())/1000, zorder=100)
ax.plot(Z.T[0], Z.T[1], zorder=1)
ax.scatter(Z[-1][0], Z[-1][1], color='xkcd:fire engine red', s=np.max(X.flatten())/1000, zorder=100)
plt.show()
'''
'''
E = myranda.EnderMan(N=1000, step=lambda: np.random.choice(np.arange(100, 10100, 100), p=VP), pidx=0)
figa, ax = plt.subplots(ncols=1, nrows=1)
ax.scatter(E.T[0], E.T[1])
plt.show()
'''
'''
figa, ax = plt.subplots(ncols=1, nrows=1)
ax = plt.axes(projection='3d')
ax.plot(xs=X.T[0], ys=X.T[1], zs=X.T[2], color='xkcd:fire engine red')
ax.plot(xs=Y.T[0], ys=Y.T[1], zs=Y.T[2], color='xkcd:blurple')
ax.plot(xs=Z.T[0], ys=Z.T[1], zs=Z.T[2], color='xkcd:tangerine')
plt.show()
'''
'''
E = myranda.EnderMan(N=1000, step=lambda: 1)
figa, ax = plt.subplots(ncols=1, nrows=1)
ax = plt.axes(projection='3d')
ax.scatter(E.T[0], E.T[1], E.T[2], color='xkcd:fire engine red')
plt.show()
'''
'''
figa, ax = plt.subplots(ncols=1, nrows=1)
ax = plt.axes(projection='3d')
for i in range(n_citt):
    last = np.zeros((n_dim, 2))
    for j in range(3):
        last[j][0] = Popolo[i].T[j][-2]
        last[j][1] = Popolo[i].T[j][-1]
    ax.plot(xs=last[0], ys=last[1], zs=last[2], color=colores[i])
plt.show()
'''
'''
n, bins, patches = plt.hist(LastVote, bins=n_dim)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
col = bin_centers - min(bin_centers)
col /= max(col)
for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', pl.cm.spring(c))
plt.show()
'''