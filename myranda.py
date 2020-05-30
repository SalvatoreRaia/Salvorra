import numpy as np
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


def Stepper(n_dim, step):
    X = np.zeros((n_dim*2, n_dim))
    for i in range(n_dim):
        X[i][i] = step
        X[i+n_dim][i] = -step
    return X


def absolute_stepper(n_dim, step):
    return step*np.identity(n_dim)


def MaxwellBoltzmann(v, m, T):  # v_picco = sqrt(2kT/m), v_rms = sqrt(3kT/m)
    kB = 1.38e-23
    return (m/(2*np.pi*kB*T))**1.5 * 4*np.pi * v**2 * np.exp(-m*v**2/(2*kB*T))


def DistributoreDiMaxwell(v_max, res, num_v, T=293.15, m=1.0):  # m in uma, T in kelvin
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

    def __init__(self, n_dim=2, spatial=True, p_inerzia=None, n_steps=1000, stepper=Stepper):  # stepper(n_dim, step)
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
        if pidx is None:
            pidx = np.random.choice(np.arange(self.n_dir))
        for i in range(1, self.n_steps):
            previous = pidx
            P = inerzia(pindex=pidx, p=self.p_inerzia, n_dir=self.n_dir)
            pidx = np.random.choice(np.arange(self.n_dir), p=P)
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


class AlleUrne(object):

    def __init__(self, n_steps, years, mega_party, party_strings, random_seed=time.time()):
        self.n_steps = n_steps
        self.years = years
        self.mega_party = mega_party
        self.party_strings = party_strings
        self.random_seed = random_seed
        np.random.seed(int(self.random_seed))

    def Movimento5Random(self, parties, lr_idx, previous, SD_ign):  # previous = (bool, name, lr_idx)
        if previous[0]:
            if bool(np.random.choice(np.array([0, 1]), p=np.array([1-self.p_inerzia, self.p_inerzia]))):
                return previous[1]
            else:
                lr_idx = np.delete(lr_idx, np.where(parties == previous[1]))
                parties = np.delete(parties, np.where(parties == previous[1]))
                return self.MazzettaDiGauss(party_str=parties, lr_idx=lr_idx, my_lr=previous[2], ign=SD_ign)
        else:
            return self.MazzettaDiGauss(party_str=parties, lr_idx=lr_idx, my_lr=previous[2], ign=SD_ign)

    @staticmethod
    def MazzettaDiGauss(party_str, lr_idx, my_lr, ign):
        lr_sbj = np.zeros(lr_idx.shape)
        for i in range(lr_idx.shape[0]):
            lr_sbj[i] = np.random.normal(loc=lr_idx[i], scale=ign, size=None)
        return party_str[np.where(np.abs(lr_sbj - my_lr) == np.min(np.abs(lr_sbj - my_lr)))]

    @staticmethod
    def MazzettaUniforme(party_str, lr_idx, my_lr, ign):
        lr_sbj = np.zeros(lr_idx.shape)
        for i in range(lr_idx.shape[0]):
            lr_sbj[i] = 2 * ign * np.random.random() + lr_idx[i] - ign
        return party_str[np.where(np.abs(lr_sbj - my_lr) == np.min(np.abs(lr_sbj - my_lr)))]

    @staticmethod
    def MazzettaIperbolica(lr_idx, my_lr, k):
        a = 0
        for i in lr_idx:
            a += (i - my_lr)**(-k)
        a = 1/a
        P = np.zeros(lr_idx.shape[0])
        for i in range(P.shape[0]):
            P[i] = a * ((lr_idx[i] - my_lr)**(-k))
        return P

    @staticmethod
    def MetricaDelta(pred, true, pop, normalized=True):
        if normalized:
            return 1 - 0.5*(np.sum(np.abs(pred - true)) / pop)
        else:
            return np.sum(np.abs(pred - true))

    def ElezioniEasy(self, p_inerzia, starting_year, sigma=1):
        # mega_party = (year, array_parties_2d) in cui ogni party è caratterizzato da un int
        # array_parties_2d = (idx, n_votes, lr_idx)
        # E è 2d : anno x n_voti
        self.p_inerzia = p_inerzia
        E = np.zeros((self.n_steps, len(self.party_strings) + 3))
        year_idx = np.where(self.years == starting_year)
        E[0] = self.mega_party[year_idx][0][1]
        for i in range(1, self.n_steps):
            # print(E[i-1])
            for n in range(E[i-1].shape[0]):
                try:
                    v = int(E[i-1][n])
                except ValueError:
                    v = 0
                for _ in range(v):
                    name_p = self.party_strings[n]
                    lr_p = self.mega_party[i-1 + year_idx[0][0]][2][n]
                    bool_p = True
                    if self.mega_party[i + year_idx[0][0]][1][n] == 0:
                        bool_p = False
                    vote = self.Movimento5Random(parties=self.party_strings[np.nonzero(self.mega_party[i + year_idx[0][0]][1])],
                                                 lr_idx=self.mega_party[i + year_idx[0][0]][2][np.nonzero(self.mega_party[i + year_idx[0][0]][1])],
                                                 previous=(bool_p, name_p, lr_p), SD_ign=sigma)
                    E[i][np.where(self.party_strings == vote)] += 1
        # print('%d/%d elezioni  |  %d/%d cittadini hanno votato'
        # % (i+1, self.n_steps, c, np.sum(mega_party[i + year_idx[0][0]][1])))
        return E

    def Sfoglio(self, mega=True, E=None, starting_year=1946):
        if mega:
            sfogli = []
            for y in range(self.n_steps):
                dizio = {}
                for n, i in enumerate(self.party_strings):
                    if self.mega_party[np.where(self.years == starting_year)[0][0] + y][1][n] != 0:
                        dizio[i] = self.mega_party[np.where(self.years == starting_year)[0][0] + y][1][n]
                sfogli.append(dizio)
        else:
            if E is None:
                return None
            sfogli = []
            for y in range(E.shape[0]):
                dizio = {}
                for n, i in enumerate(self.party_strings):
                    if E[y][n] != 0:
                        dizio[i] = E[y][n]
                sfogli.append(dizio)
        return sfogli
