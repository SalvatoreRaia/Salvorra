import numpy as np
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
import pandas as pd


def LinearRegressor(X, Y, Xerr, Yerr):
	lr = LinearRegression()
	sigma_X = Xerr / 3
	w = 1 / np.square(Yerr / 3)
	lr.fit(X, Y, w.ravel())
	for _ in range(100):
		sigma_YX = lr.coef_[0] * sigma_X
		w = 1 / (np.square(Yerr / 3) + np.square(sigma_YX))
		lr.fit(X, Y, w.ravel())
	delta = np.sum(w * np.dot(w.T, np.square(X))) - np.dot(w.T, X)**2
	merr = (np.sum(w) / delta)**0.5
	qerr = (np.dot(w.T, np.square(X)) / delta)**0.5
	return lr.coef_[0], lr.intercept_, merr.flatten(), qerr.flatten(), lr.score(X, Y)


# 0 == 1mm
df = pd.read_excel('./dati_faraday.xls', header=0, sep='\s+')
V_1 = df['1R_V_2'].values[df['Operatore_1R_2'] == 'Salvo_1']
V_2 = df['1R_V_2'].values[df['Operatore_1R_2'] == 'Salvo_2']
d_1 = df['1R_Distanza_2'].values[df['Operatore_1R_2'] == 'Salvo_1']
err_V_1 = df['Err_1R_V_2'].values[df['Operatore_1R_2'] == 'Salvo_1']
d_2 = df['1R_Distanza_2'].values[df['Operatore_1R_2'] == 'Salvo_2']
err_V_2 = df['Err_1R_V_2'].values[df['Operatore_1R_2'] == 'Salvo_2']
err_d_1 = np.ones(d_1.shape)*0.05
err_d_2 = np.ones(d_2.shape)*0.05
V_3 = df['1R_V_2'].values[df['Operatore_1R_2'] == 'Salvo_3']
d_3 = df['1R_Distanza_2'].values[df['Operatore_1R_2'] == 'Salvo_3']
err_V_3 = df['Err_1R_V_2'].values[df['Operatore_1R_2'] == 'Salvo_3']
err_d_3 = np.ones(d_3.shape)*0.05

figa, ax = plt.subplots(nrows=1, ncols=1, figsize=(10.80, 10.80))
ax.errorbar(d_1, V_1, xerr=err_d_1, yerr=err_V_1, fmt='o', markersize=1.5, ecolor='xkcd:pinkish purple', color='xkcd:neon purple', elinewidth=2, label='dati sperimentali', zorder=90)
ax.errorbar(d_2, V_2, xerr=err_d_2, yerr=err_V_2, fmt='o', markersize=1.5, ecolor='xkcd:forrest green', color='xkcd:neon green', elinewidth=2, label='dati sperimentali', zorder=90)
ax.errorbar(d_3, V_3, xerr=err_d_3, yerr=err_V_3, fmt='o', markersize=1.5, ecolor='xkcd:dark orange', color='xkcd:tangerine', elinewidth=2, label='dati sperimentali', zorder=90)

'''
figa, ax = plt.subplots(nrows=1, ncols=1, figsize=(10.80, 10.80))
ax.set_title('MISURE AL VARIARE DI $\\theta$')
ax.errorbar(angolo, f_rot, xerr=err_angolo, yerr=err_rot, fmt='o', markersize=1.5, ecolor='xkcd:pinkish purple', color='xkcd:neon purple', elinewidth=2, label='dati sperimentali', zorder=90)
ax.axhline(np.max(f_rot), color='xkcd:kelly green', label='$F_{max}$')
ax.axhline(np.max(f_rot) + err_rot[0], color='xkcd:sage', linestyle='--')
ax.axhline(np.max(f_rot) - err_rot[0], color='xkcd:sage', linestyle='--')
ax.axhline(np.min(f_rot), color='xkcd:gold', label='$F_{min}$')
ax.axhline(np.min(f_rot) + err_rot[0], color='xkcd:yellow brown', linestyle='--')
ax.axhline(np.min(f_rot) - err_rot[0], color='xkcd:yellow brown', linestyle='--')
ax.set_xticks(np.arange(-120, 120, 1), minor=True)
ax.set_xlim((-120, 120))
ax.set_yticks(np.arange(104.5, 105.6, 0.01), minor=True)
ax.set_xticks(np.arange(-120, 125, 15), minor=False)
ax.set_yticks(np.arange(104.5, 105.65, 0.1), minor=False)
ax.grid(color='xkcd:light violet', linestyle=':', which='minor', alpha=0.75)
ax.grid(color='xkcd:orchid', linestyle='--', which='major', zorder=1, alpha=0.75)
ax.set_xlabel('$\\theta$ (°)')
ax.set_ylabel('Forza totale ($g$)')
ax.legend(loc='lower right')
plt.savefig('./EM_grafico_4.svg')
'''
