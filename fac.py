import numpy as np
import pandas as pd
from itertools import groupby
from collections import Counter
import matplotlib.pyplot as plt
import pickle


#IMPORTING DATA
with open(r'C:\Users\Salvo\Desktop\y.pkl', 'rb') as f:
    y_tot = pickle.load(f)

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)

#total no. of counties
N_counties = len(Adata[1][:]) -1

#choose the county by a number between [1,N_counties]
#c = np.random.choice(np.arange(1, N_counties + 1 ))
#print(c)
c = 1712
c1 = 273
c2 = 2585
county = Adata[1][c]
county1 = Adata[1][c1]
county2 = Adata[1][c2]


#chooose the number of words you need; default is the maximum
word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False, index_col=[0])
word0 = int(word.iloc[c, 1])
word1 = int(word.iloc[c1, 1])
word2 = int(word.iloc[c2, 1])

print(word, word0, word1, word2)

#IMPORTING DATA
filename = r"C:\Users\Salvo\Desktop\IFISC\data\r_f_county_w\freq_ranking_"+county+'(' + str(c) +')_'+str(word0)+"w.csv"
filename1 = r"C:\Users\Salvo\Desktop\IFISC\data\r_f_county_w\freq_ranking_"+county1+'(' + str(c1) +')_'+str(word1)+"w.csv"
filename2 = r"C:\Users\Salvo\Desktop\IFISC\data\r_f_county_w\freq_ranking_"+county2+'(' + str(c2) +')_'+str(word2)+"w.csv"

dat = pd.read_csv(filename)
dat1 = pd.read_csv(filename1)
dat2 = pd.read_csv(filename2)

print(dat)
print(dat1)
print(dat2)


dat = dat.sort_values(by=['words'], ascending=True)
dat1 = dat.sort_values(by=['words'], ascending=True)
dat2 = dat.sort_values(by=['words'], ascending=True)

print(dat)
print(dat1)
print(dat2)

w = dat.iloc[:, 1].values.tolist()
f = dat.iloc[:, 2].values.tolist()

w1 = dat1.iloc[:, 1].values.tolist()
f1 = dat1.iloc[:, 2].values.tolist()

w2 = dat2.iloc[:, 1].values.tolist()
f2 = dat2.iloc[:, 2].values.tolist()

#print(w, f)

w = w + w1 + w2
f = f + f1 + f2
print(len(w), len(f))

#wordcount = [len(list(group)) for key, group in groupby(w)]

#x = 3
#d = Counter(l)
wordcount2 = Counter(w)



mergedcounties = pd.DataFrame({'words': w, 'frequency': f})
mergedcounties= mergedcounties.sort_values(by=['words'], ascending=True)
print(mergedcounties)





f = mergedcounties.iloc[:, 1].values.tolist()
print(len(f))

f_merged = []

myset = set(w)
#print(myset)
uniquewords = list(myset)

uniquewords = pd.DataFrame({'words': uniquewords})
uniquewords= uniquewords.sort_values(by=['words'], ascending=True)
uniquewords = uniquewords.iloc[:, 0].values.tolist()


for i in range(len(uniquewords)):
    templist = [np.sum(f[i:i+wordcount2[uniquewords[i]]])]
    f_merged += templist
    #print(f_merged)

print(f_merged)
print(uniquewords)
print(uniquewords[0], uniquewords[len(uniquewords)-1])


f_merged = pd.DataFrame({'f': f_merged})
f_merged= f_merged.sort_values(by=['f'], ascending=False)
f_merged = f_merged.iloc[:, 0].values.tolist()

print(f_merged)

y = np.int64(f_merged)
sum_y = np.sum(y)
#Normalizing
y = y/sum_y
norm_sum_y = np.sum(y)
print('Are the data normalized ?:', norm_sum_y)
print(y_tot[1])

x = np.arange(1, len(f_merged)+1)



# GRAPHIC
fig = plt.figure(figsize=[2.,3.])
ax1 = fig.add_subplot(1, 1, 1)

ax1.scatter(x, y,
            s=0.1, label='The 3 merged counties', color='xkcd:royal blue')

ax1.scatter(np.arange(1, len(y_tot[c-1])+1), y_tot[c-1],
            s=0.1, label=str(c), color = 'red')

ax1.scatter(np.arange(1, len(y_tot[c1-1])+1), y_tot[c1-1],
            s=0.1, label=str(c1))

ax1.scatter(np.arange(1, len(y_tot[c2-1])+1), y_tot[c2-1],
            s=0.1, label=str(c2))
x1=np.arange(2, 20000)
ax1.plot(x1, -np.log10(x1)**2)

#Axis grid
ax1.grid(color='xkcd:light blue', linestyle='--', which='major', zorder=1, alpha=0.85)
ax1.grid(color='xkcd:light blue', linestyle='--', which='minor', zorder=1, alpha=0.35)

ax1.set_xticks(np.arange(1,100, 10), minor=True)


# Figure settings

# Setting log-log scale
ax1.set_yscale('log')
ax1.set_xscale('log')

# Axis titles
ax1.set_xlabel('Ranking')
ax1.set_ylabel('Normalized frequency')

ax1.legend(loc='best')
ax1.set_ylim([10**-7,1])


plt.show()



