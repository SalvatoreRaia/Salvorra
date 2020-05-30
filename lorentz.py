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
	qerr = (np.dot(w.T, np.square(X))/ delta)**0.5
	return lr.coef_[0], lr.intercept_, merr.flatten(), qerr.flatten(), lr.score(X, Y)


df = pd.read_excel('./dati_lorentz.xls', header=0, sep='\s+')
F = df[['Forza peso (g)']].values[:15]
I = df[['Corrente A']].values[:15]
Ferr = df[['errore forza']].values[:15]
Ierr = df[['errore I']].values[:15]
F_2 = df[['Forza peso (g)']].values[15:32]
I_2 = df[['Corrente A']].values[15:32]
Ferr_2 = df[['errore forza']].values[15:32]
Ierr_2 = df[['errore I']].values[15:32]
F_3 = df[['Forza peso (g)']].values[32:37]
I_3 = df[['Corrente A']].values[32:37]
Ferr_3 = df[['errore forza']].values[32:37]
Ierr_3 = df[['errore I']].values[32:37]

f_rot = df[['Forza peso rot (g)']].values
angolo = df[['Angolo (deg)']].values
err_rot = df[['errore forza 3']].values
err_angolo = df[['errore angolo']].values

out_len = df[['vari fili out len']][:6].values
int_len = df[['vari fili int len']][:6].values
len_fili = 0.5*(out_len + int_len).flatten()
delta_f = df[['vari fili 2A']][:6].values.flatten() - df[['vari fili 0A']][:6].values.flatten()
err_out = np.array([np.max(df[['lunghezza (cm) esterna 1 39']][:5].values) - np.min(df[['lunghezza (cm) esterna 1 39']][:5].values),
					np.max(df[['lunghezza (cm) esterna 2 37']][:3].values) - np.min(df[['lunghezza (cm) esterna 2 37']][:3].values),
					np.max(df[['lunghezza (cm) esterna 3 38']][:3].values) - np.min(df[['lunghezza (cm) esterna 3 38']][:3].values),
					np.max(df[['lunghezza (cm) esterna 4 40']][:3].values) - np.min(df[['lunghezza (cm) esterna 4 40']][:3].values),
					np.max(df[['lunghezza (cm) esterna 5 41 (X2)']][:3].values) - np.min(df[['lunghezza (cm) esterna 5 41 (X2)']][:3].values),
					np.max(df[['lunghezza (cm) esterna 6 42 (x2)']][:3].values) - np.min(df[['lunghezza (cm) esterna 6 42 (x2)']][:3].values)])/2.
err_int = np.array([np.max(df[['lunghezza (cm) interna 1 39']][:5].values) - np.min(df[['lunghezza (cm) interna 1 39']][:5].values),
					np.max(df[['lunghezza (cm) interna 2 37']][:3].values) - np.min(df[['lunghezza (cm) interna 2 37']][:3].values),
					np.max(df[['lunghezza (cm) interna 3 38']][:3].values) - np.min(df[['lunghezza (cm) interna 3 38']][:3].values),
					np.max(df[['lunghezza (cm) interna 4 40']][:3].values) - np.min(df[['lunghezza (cm) interna 4 40']][:3].values),
					np.max(df[['lunghezza (cm) interna 5 41 (X2)']][:3].values) - np.min(df[['lunghezza (cm) interna 5 41 (X2)']][:3].values),
					np.max(df[['lunghezza (cm) interna 6 42 (x2)']][:3].values) - np.min(df[['lunghezza (cm) interna 6 42 (x2)']][:3].values)])/2.
err_fili = (err_out + err_int)*0.5
err_delta_f = df[['errore forza']].values[:6]
m_f, q_f, m_err_f, q_err_f, lincorr_f = LinearRegressor(len_fili.reshape(-1, 1), delta_f.reshape(-1, 1), err_fili.reshape(-1, 1), err_delta_f)

tara_i = df[['peso ideale']][:8].values
tara_r = df[['peso reale']][:8].values
err_tr = df[['errore reale']][:8].values
err_ti = df[['errore peso ideale']][:8].values

m, q, m_err, q_err, lincorr = LinearRegressor(I, F, Ierr, Ferr)# gli err sono sigma!!!
m_2, q_2, m_err_2, q_err_2, lincorr_2 = LinearRegressor(I_2, F_2, Ierr_2, Ferr_2)
m_3, q_3, m_err_3, q_err_3, lincorr_3 = LinearRegressor(I_3, F_3, Ierr_3, Ferr_3)
m_err, q_err = m_err*3, q_err*3
m_err_2, q_err_2 = m_err_2*3, q_err_2*3
m_err_3, q_err_3 = m_err_3*3, q_err_3*3
print('\nCoefficienti set 1:\n', m, q, m_err, q_err, lincorr)
print('\n\nCoefficienti set 2:\n',m_2, q_2, m_err_2, q_err_2, lincorr_2)
print('\n\nCoefficienti set 3:\n',m_3, q_3, m_err_3, q_err_3, lincorr_3)
print('\n\nCoefficienti varie len:\n',m_f, q_f, m_err_f, q_err_f, lincorr_f)
print('\n\nMassimo sinusoide:\t', np.max(f_rot), err_rot[0])
print('\nMinimo sinusoide:\t', np.min(f_rot), err_rot[0])

'''
X = np.arange(-0.1, 5.6, 0.5)
figa, ax = plt.subplots(nrows=1, ncols=2, figsize=(19.20, 10.80))
ax[0].set_title('PRIMO SET DI MISURE (r = 0.9885)')
ax[0].errorbar(I, F, xerr=Ierr, yerr=Ferr, fmt='o', markersize=1.5, ecolor='xkcd:pinkish red', color='xkcd:scarlet', elinewidth=2, label='dati sperimentali', zorder=90)
ax[0].plot(X, m*X + q, label='$F_{tot}$ = -0.085*I + 155.652', zorder=100)# label da modificare coi numeri
ax[0].fill_between(X, (m-m_err)*X + q-q_err, (m+m_err)*X + q+q_err, color='xkcd:pale blue', label="zona nell'errore della retta")
ax[0].set_xticks(np.arange(-0.1, 5.1, 0.1), minor=True)
ax[0].set_yticks(np.arange(155.2, 155.75, 0.025), minor=True)
ax[0].set_xticks(np.arange(-0.1, 5.5, 0.5), minor=False)
ax[0].set_yticks(np.arange(155.2, 155.75, 0.1), minor=False)
ax[0].set_xlim((-0.25, 5.2))
ax[0].set_ylim((155.15, 155.8))
ax[0].set_xlabel('Intensità di corrente ($A$)')
ax[0].set_ylabel('Forza totale ($g$)')
ax[0].grid(color='xkcd:baby blue', linestyle=':', which='minor', alpha=0.75)
ax[0].grid(color='xkcd:sea blue', linestyle='--', which='major', zorder=1, alpha=0.75)
ax[0].legend(loc='upper right')

ax[1].set_title('SECONDO SET DI MISURE (r = 0.9859)')
ax[1].errorbar(I_2, F_2, xerr=Ierr_2, yerr=Ferr_2, fmt='o', markersize=1.5, ecolor='xkcd:pinkish red', color='xkcd:scarlet',elinewidth=2, label='dati sperimentali', zorder=90)
# ax[1].errorbar(I_3, F_3, xerr=Ierr_3, yerr=Ferr_3, fmt='o', markersize=1.5, ecolor='xkcd:blurple', color='xkcd:scarlet', elinewidth=2, label='dati sperimentali', zorder=90)
ax[1].plot(X, m_2*X + q_2, label='$F_{tot}$ = -0.073*I + 155.663', zorder=100, color='xkcd:leaf green')
ax[1].fill_between(X, (m_2-m_err_2)*X + q_2-q_err_2, (m_2+m_err_2)*X + q_2+q_err_2, color='xkcd:very light green', label="zona nell'errore della retta", alpha=0.5)
ax[1].set_xticks(np.arange(-0.1, 5.1, 0.1), minor=True)
ax[1].set_yticks(np.arange(155.2, 155.75, 0.025), minor=True)
ax[1].set_xticks(np.arange(-0.1, 5.5, 0.5), minor=False)
ax[1].set_yticks(np.arange(155.2, 155.75, 0.1), minor=False)
ax[1].set_xlim((-0.25, 5.2))
ax[1].set_ylim((155.15, 155.8))
ax[1].set_xlabel('Intensità di corrente ($A$)')
ax[1].set_ylabel('Forza totale ($g$)')
ax[1].grid(color='xkcd:chartreuse', linestyle=':', which='minor')
ax[1].grid(color='xkcd:jungle green', linestyle='--', which='major', zorder=1, alpha=0.75)
ax[1].legend(loc='upper right')
plt.savefig('./EM_grafico_1.svg')


X = np.arange(-0.1, 5.6, 0.5)
figa, ax = plt.subplots(nrows=1, ncols=1, figsize=(10.80, 10.80))
ax.set_title('CONFRONTO TRA I SET')
ax.errorbar(I, F, xerr=Ierr, yerr=Ferr, fmt='o', markersize=1.5, ecolor='xkcd:sea blue', color='xkcd:scarlet', elinewidth=2, label='primo set', zorder=90)
ax.errorbar(I_2, F_2, xerr=Ierr_2, yerr=Ferr_2, fmt='o', markersize=1.5, ecolor='xkcd:leaf green', color='xkcd:scarlet', elinewidth=2, label='secondo set', zorder=90)
ax.errorbar(I_3, F_3, xerr=Ierr_3, yerr=Ferr_3, fmt='o', markersize=1.5, ecolor='xkcd:pinkish red', color='xkcd:scarlet', elinewidth=2, label='set decisivo (r = 0.9981)', zorder=90)
ax.set_xticks(np.arange(-0.1, 5.1, 0.1), minor=True)
ax.set_yticks(np.arange(155.2, 155.75, 0.025), minor=True)
ax.set_xticks(np.arange(-0.1, 5.5, 0.5), minor=False)
ax.set_yticks(np.arange(155.2, 155.75, 0.1), minor=False)
ax.set_xlim((-0.25, 5.2))
ax.set_ylim((155.15, 155.8))
ax.set_xlabel('Intensità di corrente ($A$)')
ax.set_ylabel('Forza totale ($g$)')
ax.grid(color='xkcd:maize', linestyle=':', which='minor')
ax.grid(color='xkcd:burnt yellow', linestyle='--', which='major', zorder=1, alpha=0.75)
ax.legend(loc='upper right')
plt.savefig('./EM_grafico_2.svg')

X = np.arange(-0.1, 8.6, 0.5)
figa, ax = plt.subplots(nrows=1, ncols=1, figsize=(10.80, 10.80))
ax.set_title('MISURE AL VARIARE DI $\mathscr{l}$  (r = 0.9990)')
ax.errorbar(len_fili, delta_f, xerr=err_fili, yerr=err_delta_f, fmt='o', markersize=1.5, ecolor='xkcd:pinkish red', color='xkcd:scarlet', elinewidth=2, label='dati sperimentali', zorder=90)
ax.plot(X, m_f*X + q_f, label='$F_{L}$ = -0.0516*I + 0.017', zorder=100, color='xkcd:orange')
ax.fill_between(X, (m_f-m_err_f)*X + q_f-q_err_f, (m_f+m_err_f)*X + q_f+q_err_f, color='xkcd:tangerine', label="zona nell'errore della retta", alpha=0.15)
ax.set_xticks(np.arange(-0.1, 8.6, 0.1), minor=True)
ax.set_yticks(np.arange(-0.5, 0.1, 0.01), minor=True)
ax.set_xticks(np.arange(-0.1, 8.6, 0.5), minor=False)
ax.set_yticks(np.arange(-0.5, 0.1, 0.05), minor=False)
ax.grid(color='xkcd:peach', linestyle=':', which='minor')
ax.grid(color='xkcd:dark orange', linestyle='--', which='major', zorder=1, alpha=0.75)
ax.set_xlabel('Intensità di corrente ($A$)')
ax.set_ylabel('Forza di Lorentz ($g$)')
ax.legend(loc='upper right')
plt.savefig('./EM_grafico_3.svg')
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