import numpy as np
import pandas as pd
import time
import myranda as myr
import matplotlib.pyplot as plt

filez = np.load('./mp_pn_ys_ita_easy_10000.npz')
Italia = myr.AlleUrne(n_steps=4, years=filez['arr_2'], mega_party=filez['arr_0'], party_strings=filez['arr_1'], random_seed=8765689)
voto = Italia.ElezioniEasy(p_inerzia=0.9, starting_year=2006, R_ign=1.5)
# voto = Italia.Sfoglio(mega=False, E=voto)
voto_vero = Italia.Sfoglio(starting_year=2006)
DV = np.zeros(voto.shape)
for i, v in enumerate(voto):
    DV[i] = np.abs(v - Italia.mega_party[np.where(Italia.years == 2006)[0][0] + i][1])*0.5 / np.sum(Italia.mega_party[i + np.where(Italia.years == 2006)[0][0]][1])

c = ['xkcd:fire engine red', 'xkcd:snot green', 'xkcd:ocean blue']
figa, ax = plt.subplots(nrows=3, ncols=1)
for g in range(1, voto.shape[0]):
    ax[g-1].barh(y=filez['arr_1'][np.where(DV[g] != 0)], width=DV[g][np.where(DV[g] != 0)], color=c[g-1])
plt.show()