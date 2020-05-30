import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import powerlaw


def sigma(alpha, n):
    """
    Clauset et al 2007 equation 3.2:
        sigma = (alpha-1)/sqrt(n)
    """
    return (alpha-1.) / n**0.5


def discrete_alpha_mle(data, xmin):
    """
    Equation B.17 of Clauset et al 2009

    The Maximum Likelihood Estimator of the "scaling parameter" alpha in the
    discrete case is similar to that in the continuous case
    """
    # boolean indices of positive data
    gexmin = (data>=xmin)
    nn = gexmin.sum()
    if nn < 2:
        return 0
    xx = data[gexmin]
    alpha = 1.0 + float(nn) * (sum(np.log(xx/(float(xmin)-0.5))))**-1
    return alpha


def discrete_ksD(data, xmin, alpha):
    """
    given a sorted data set, a minimum, and an alpha, returns the power law ks-test
    D value w/data

    The returned value is the "D" parameter in the ks test

    (this is implemented differently from the continuous version because there
    are potentially multiple identical points that need comparison to the power
    law)
    """

    zz = np.sort(data[data >= xmin])
    nn = float(len(zz))
    if nn < 2:
        return np.inf
    # cx = np.arange(nn,dtype='float')/float(nn)
    # cf = 1.0-(zz/xmin)**(1.0-alpha)
    model_cdf = 1.0 - (zz.astype('float') / float(xmin)) ** (1.0 - alpha)
    data_cdf = np.searchsorted(zz, zz, side='left') / (float(nn))
    ks = max(abs(data_cdf - model_cdf))
    return ks

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)

#total no. of counties
N_counties = len(Adata[1][:]) -1

#choose the county by a number between [1,N_counties]
#c = np.random.choice(np.arange(1, N_counties + 1 ))
#print(c)
c = 2965 #custom choice
county = Adata[1][c]

#chooose the number of words you need; default is the maximum
word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False, index_col=[0])
word = int(word.iloc[c, 1])
print(word)

#IMPORTING DATA
filename = r"C:\Users\Salvo\Desktop\IFISC\data\r_f_county_w\freq_ranking_"+county+'(' + str(c) +')_'+str(word)+"w.csv"
dat = pd.read_csv(filename)
print(dat)


b = word

#Real normalized data
#PMF
x = np.arange(1, b + 1)
y = np.int64(dat.iloc[:, 2])
sum_y = np.sum(y)
print('sum y',sum_y)
y = y/sum_y
norm_sum_y = np.sum(y)
print('Are the data normalized ?:', norm_sum_y)


print('NON binned data :')
#GETTING NON BINNED DATA
y = np.int64(dat.iloc[:, 2])
dat_no_binned = np.ones(int(10**8/2), dtype= np.int8).tolist()

#array of the No. of times a word is said (each word)
tokens_sum = np.sum(dat.iloc[:,2])
w_times = np.zeros(word)
tot_w_times = 0
'''
for i in range (word):
    w_times[i] = dat.iloc[i,2] * 10**9 / tokens_sum
    #print(w_times[i], word -i)
    tot_w_times = tot_w_times + w_times[i]

for i in range (len(w_times)):
    #print(len(w_times)-i)
    values = (i+1)*np.ones(int(w_times[i]/(w_times[len(w_times)-1]))) +1
    dat_no_binned.append(values)
    print(len(w_times)-i)

print(type(dat_no_binned))
'''
#np.array(dat_no_binned, dtype=np.int8)



#Analisys on non binned data
fit = powerlaw.Fit(dat_no_binned, discrete = True)
alpha = fit.power_law.alpha
print('powerlaw fit, alpha = ', alpha)
xmin = fit.xmin
print('powerlaw fit, xmin = ', xmin)


#Survival function: sto trovando prima cdf. cdf[i] è la somma di tutti meno la somma da i fino a len(cdf)
# rimane la somma da 0 a i in posizione i. Poi faccio 1 -cdf
sf_y = np.ones(word)
for i in range (word-1):
    sf_y[i+1] = (sum_y - np.sum(y[i+1:word]))/sum_y
    #print(sf_y[i], i)

#sf_y[word] =  norm_sum_y
#sf_y = 1 - sf_y


#GRAPHIC
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.set_title('asdfghl')

ax1.plot(np.arange(len(y)),y, linewidth = 0.5, label=str(c)+' '+county+', '+str(word)+' words.', color = 'red')
#ax1.plot(np.arange(len(sf_y)),sf_y, linewidth = 0.5, label='CDF')
#ax1.plot(np.arange(len(y)),y, linewidth = 0.5, color = 'red')


#Setting log-log scale
ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.set_xlabel('fhgcjhj')
ax1.set_ylabel('dfgchv')
#ax1.grid(color='xkcd:cherry', linestyle='--', which='major', zorder=1, alpha=0.75)
#ax1.set_xticks(np.arange(0, word + 1), minor=False)



ax1.legend(loc='best')
plt.show()