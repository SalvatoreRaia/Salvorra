import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pylab as pl
import time

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

    def Movimento5Random(self, parties, lr_idx, previous, mazz_k):# previous = (bool, name, lr_idx)
        if previous[0]:
            if bool(np.random.choice(np.array([0, 1]), p=np.array([1-self.p_inerzia, self.p_inerzia]))):
                return party_space[np.where(parties == previous[1])]
            else:
                np.delete(lr_idx, np.where(parties == previous[1]))
                np.delete(party_space, np.where(parties == previous[1]))
                return np.random.choice(parties, p=self.Mazzetta(lr_idx=lr_idx, my_lr=previous[2], k=mazz_k))
        else:
            return np.random.choice(parties, p=self.Mazzetta(lr_idx=lr_idx, my_lr=previous[2], k=mazz_k))

    def Mazzetta(self, lr_idx, my_lr, k):
        a = 0
        for i in lr_idx:
            a += (i - my_lr)**(-k)
        a = 1/a
        P = np.zeros(lr_idx.shape[0])
        for i in range(P.shape[0]):
            P[i] = a * ((lr_idx[i] - my_lr)**(-k))
        return P

    def Elezioni(self, n_citt, mega_party, party_strings, E_0, random_seed=time.time()):
        # mega_party = (year, array_parties_2d) in cui ogni party è caratterizzato da un int
        # array_parties_2d = (idx, n_votes, lr_idx)
        # E è 2d : anno x n_voti
        np.random.seed(int(random_seed))
        E = np.zeros((n_steps, len(party_strings)))
        E[0] = E_0
        for i in range(1, self.n_steps):
            for n in range(E[i-1].shape):
                while E[i-1][n] != 0:
                    name_p = party_strings[n]
                    lr_p = mega_party[i-1][2]
                    bool_p = True
                    if E[i][n] == 0:
                        bool_p = False
                    vote = self.Movimento5Random(parties=party_strings[np.nonzero(mega_party[i][1])],
                                                 lr_idx=mega_party[i][2][np.nonzero(mega_party[i][1])],
                                                 previous=(bool_p, name_p, lr_p), mazz_k=1)
                    E[i][party_strings.index(vote)] += 1
                    E[i-1][n] -= 1
                print('%d/%d elezioni  |  %d/%d cittadini hanno votato' % (i, self.n_step, n, self.n_citt))
        print('RUSPAAAAA!!!')
        return E


# VP = DistributoreDiMaxwell(10000, 0.01, 100)
n_citt = 1000
n_steps = 20
n_dim = 42
myranda = RandomManhattan(n_dim=n_dim, spatial=False, p_inerzia=0.25, n_steps=n_steps, stepper=stepper)

Popolo = np.zeros((n_citt, n_steps, n_dim))
LastVote = np.zeros(n_citt)
for n in range(Popolo.shape[0]):
    Popolo[n], LastVote[n] = myranda.LukeRandomwalker(step=lambda: 1, random_seed=time.time()+n, pidx=39)
    print('Faccio propaganda: ', n+1, '/', n_citt)
colores = pl.cm.gist_rainbow(np.linspace(0, 1, n_citt))

# X = myranda.LukeRandomwalker(step=lambda: 1, pidx=0)
# Y = myranda.LukeRandomwalker(step=lambda: 1, random_seed=3450954, pidx=1)
# Z = myranda.LukeRandomwalker(step=lambda: 1, random_seed=2057, pidx=2)
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

n, bins, patches = plt.hist(LastVote, bins=n_dim)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
col = bin_centers - min(bin_centers)
col /= max(col)
for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', pl.cm.spring(c))
plt.show()
