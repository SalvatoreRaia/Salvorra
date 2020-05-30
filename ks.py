import numpy as np
import pandas as pd
import pickle

Adata = pd.read_csv('/home/salvatore/Project1/INDEXES_PBW/INDEX_A_PBW.txt', sep=",", header=None, low_memory=False)

# total no. of counties
TOT_counties = len(Adata[1][:]) - 1

# List of total number of words per county
tot_word = pd.read_csv('/home/salvatore/Project1/no_words_per_county.csv', sep=",", low_memory=False, index_col=[0])
print(tot_word)

c = np.arange(0, TOT_counties)
word = np.zeros(TOT_counties)

print('ciao')



#Loadind y
which_county = np.arange(1, TOT_counties )
with open('/home/salvatore/Project1/y.pkl', 'rb') as f:
    y = pickle.load(f)

print('ciao ora inizia il for I')
cdf_y = [[]]
# CDF
for i in range(1, TOT_counties):
    word[i] = int(tot_word.iloc[c[i], 1])
    b = int(word[i])
    cdf_y.append([])
    partial_sum = 0
    # print(b, len(y[i]))

    for j in range(b):
        partial_sum = partial_sum + y[i - 1][j]
        cdf_y[i].append(partial_sum)

print('CDF calculated')

print('a check', cdf_y[0], cdf_y[1])

# KS Statistics
grid = np.arange(1, TOT_counties)
myD = np.zeros((len(grid), len(grid)))

# Filling triangular matrix
"""
Probably the filling doesn't include the last county.
I added it manually afterward. This should be fixed but the distance calulation for
a single couple of counties is ok
"""
for i in range(len(grid)):
    for j in range(i, len(grid)):
        minNWords = int(min(word[i + 1], word[j + 1]))
        myD[i][j] = np.amax(abs(np.array(cdf_y[i + 1][:minNWords]) - np.array(cdf_y[j + 1][:minNWords])))

    print(len(grid) - i)

combinations = len(grid) * (len(grid) - 1) * 0.5
print(combinations)
np.savez('/home/salvatore/Project1/all_counties/ks/run/myD_value_' + str(int(combinations)) + '_combinations', D=myD,
         counties=grid)
print('Done')





