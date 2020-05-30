import numpy as np
import pandas as pd
from scipy.stats import zipf
from scipy import stats

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)

#total no. of counties
N_counties = len(Adata[1][:]) -1
print(N_counties)

#choose the county by a number between [1,N_counties], (random choice or user setting)
c = np.random.choice(np.arange(1, N_counties + 1 ))
print(c)
c = 1968
county = Adata[1][c]

#chooose the number of words you need; default is the maximum
word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False, index_col=[0])
word = int(word.iloc[c, 1])
print(word)

#IMPORTING DATA
filename = r"C:\Users\Salvo\Desktop\IFISC\data\r_f_county_w\freq_ranking_"+county+'('+str(c)+')_'+str(word)+"w.csv"
dat = pd.read_csv(filename)
print('dat',dat)

#array of the No. of times a word is said (each word)
tokens_sum = np.sum(dat.iloc[:,2])
print('tockens sum',tokens_sum)
w_times = np.zeros(word)
tot_w_times = 0
for i in range (word):
    w_times[i] = dat.iloc[i,2] * 10**9 / tokens_sum
    tot_w_times = tot_w_times + w_times[i]
    print('w_times[]',w_times[i])

print('w_times',w_times)
print('sumw_times', np.sum(w_times)/10**9)
print('tot_w_times',tot_w_times)
print('Check tot_w', np.sum(tot_w_times))

rel_freq = w_times/tot_w_times
print(rel_freq)
a = np.sum(rel_freq)
print(a)
b = np.sum(dat.iloc[:,2]/tokens_sum)
print(b)
