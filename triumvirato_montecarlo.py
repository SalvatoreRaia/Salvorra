import numpy as np
import pandas as pd
import time
import myranda as myr
import matplotlib.pyplot as plt
import seaborn as sbor

nation = 'spa'
filez = np.load('./mp_pn_ys_' + nation + '_easy_10000.npz')
yearray = filez['arr_2']
Italia = myr.AlleUrne(n_steps=2, years=filez['arr_2'], mega_party=filez['arr_0'], party_strings=filez['arr_1'],
                      random_seed=time.time())
very_beginning = time.time()
PB = np.zeros(yearray.shape[0])
SB = np.zeros(yearray.shape[0])
mcn = (256, 64, 16, 4)
p_ray = (1, 0.1, 0.001, 0.0001)
s_ray = (5, 0.5, 0.005, 0.0005)
runs = 14
done = 0

# Gli intervalli sono: Paride_1 = [:7], Paride_2 = [6:-7], Salvo = [-8:] da mettere al posto delle [:]
# Inoltre crea una cartella che si chiama così:  logs_mc_p_sigma
truncated_yearray = (np.arange(yearray.shape[0])[8:20], yearray[8:20])

print('Stai iniziando il processo con ' + str(truncated_yearray[0].shape[0]) + ' anni (intervallo ' +
      str(truncated_yearray[1][0]) + '-' + str(truncated_yearray[1][-1]) + ')\n')

for y in range(truncated_yearray[0].shape[0] - 1):
    data = 'Anno ' + str(truncated_yearray[1][y]) + ' | '
    log_path = './logs_mc_p_sigma/' + nation + '/anno' + str(truncated_yearray[1][y]) + '.txt'
    fp = open(log_path, 'w+')
    p = None
    s = None
    p_range = (0.0001, 0.9999)
    s_range = (0, 5)
    for i in range(1, len(p_ray) + 1):
        D = np.zeros(mcn[i-1])
        P = np.zeros(mcn[i-1])
        S = np.zeros(mcn[i-1])
        iters = 'iter ' + str(i) + '/4 | '
        for j in range(mcn[i-1]):
            P[j] = (p_range[1] - p_range[0]) * np.random.random_sample() + p_range[0]
            S[j] = (s_range[1] - s_range[0]) * np.random.random_sample() + s_range[0]
            D[j] = np.mean(np.array([Italia.MetricaDelta(
                Italia.ElezioniEasy(p_inerzia=P[j], starting_year=truncated_yearray[1][y], sigma=S[j])[1],
                filez['arr_0'][truncated_yearray[0][y+1]][1], pop=10000) for _ in range(runs)]))
            done += 1
            print(data + iters + 'montecarlo ' + str(j+1) + '/' + str(mcn[i-1]) +
                  ' | p = ' + str(p) + ' | sigma = ' + str(s) + '\nValori temporanei: D = ' + str(D[j]) + ' | p = ' +
                  str(P[j]) + ' | sigma = ' + str(S[j]))
            fp.write(data + iters + 'montecarlo ' + str(j+1) + '/' + str(mcn[i-1]) + ' | p = ' + str(p) +
                     ' | sigma = ' + str(s) + '\nValori temporanei: D = ' + str(D[j]) + ' | p = ' + str(P[j]) +
                     ' | sigma = ' + str(S[j]) + '\n')
        d = np.max(D)
        p = P[D == d][0]
        s = S[D == d][0]
        try:
            if p - p_ray[i] > 0 and p + p_ray[i] < 1:
                p_range = (p - p_ray[i], p + p_ray[i])
            elif p - p_ray[i] <= 0 and p + p_ray[i] < 1:
                p_range = (0.0001, p + p_ray[i])
            elif p - p_ray[i] > 0 and p + p_ray[i] >= 1:
                p_range = (p - p_ray[i], 0.9999)
            else:
                p_range = (0.0001, 0.9999)
            if s - s_ray[i] < 0:
                s_range = (0, s + s_ray[i])
            else:
                s_range = (s - s_ray[i], s + s_ray[i])
        except IndexError:
            pass
    PB[y] = p
    SB[y] = s
    fp.close()
np.savez('./valori_mc_p_sigma/valori_mc_p_sigma_' + nation + '_' + str(truncated_yearray[1][0]) + '-' + str(truncated_yearray[1][-1]) + '.npz', PB, SB)
print('\nDopo poco più di ' + str((time.time() - very_beginning)//3600) + ' ore hai ottenuto i seguenti valori:')
print(PB)
print(SB)
