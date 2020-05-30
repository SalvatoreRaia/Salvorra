import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zipf

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)

#total no. of counties
N_counties = len(Adata[1][:]) -1

#choose the county by a number between [1,N_counties]
c = np.random.choice(np.arange(1, N_counties + 1 ))
c = 2345
county = Adata[1][c]

#chooose the number of words you need; default is the maximum
word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False, index_col=[0])
word = int(word.iloc[c, 1])
print(word)

#IMPORTING DATA
filename = r"C:\Users\Salvo\Desktop\IFISC\data\r_f_county_w\freq_ranking_"+county+'(' + str(c) +')_'+str(word)+"w.csv"
dat = pd.read_csv(filename)
print(dat)

alpha = 1.5
random_zipf = zipf.pmf(zipf.rvs(a = alpha, size=1000), a = alpha)
random_zipf = -np.sort(-random_zipf)
print(random_zipf)


# GRAPHIC
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)


x = np.arange(1,int(len(random_zipf))+1)
print(len(x), len(random_zipf))
#ax1.plot(x, y, label=county, color='red', linewidth = 1)
ax1.scatter(x,random_zipf, label="Zipf's sampled data", s=1, color = 'red')
#ax1.scatter(x_sample,random_zipf1,  label="Zipf's scipy sampled data", s=1, color='orange')

ax1.grid(color='xkcd:cherry', linestyle='--', which='major', zorder=1, alpha=0.75)
ax1.set_xticks(np.arange(0, word + 1), minor=False)

#Setting log-log scale
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend(loc='upper right')

plt.show()