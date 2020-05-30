import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import scipy.stats as stats

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)

#total no. of counties
TOT_counties = len(Adata[1][:]) -1

# List of total number of words per county
tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,index_col=[0])
print(tot_word)
# How many counties you want to analyse? Change first number
n_counties = 12 + 1
S = 1

#number of words of the i-th county
word = np.zeros(n_counties)
#indx o i-th county
c = np.arange(n_counties)
print('ciao')
y = [0]
#cdf_y = []
#Initializing variables
x = np.arange(1, 70000)
for i in range (S,n_counties):
    word[i] = int(tot_word.iloc[c[i], 1])
    b = int(word[i])

    y.append([np.zeros(b)])


cdf_y = y
print(y)
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

for i in range (S,n_counties):
    #word[i] = int(tot_word.iloc[c[i], 1])
    print(word[i])
    print(Adata[1][c[i]])
    print(y[i])

    #IMPORTING DATA
    filename = "C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\\row_frequencies\\frequencies_" + str(Adata[1][i]) + "_(" + str(c[i]) + ")" + ".txt"
    dat = pd.read_csv(filename, sep=" ", header=None)
    #print(dat)

    b = int(word[i])

    #Real normalized data
    #PMF
    y[i] = dat.iloc[:,0]
    y[i] = - np.sort(-y[i])
    y[i] = y[i]/np.sum(y[i])
    print('Are the data normalized ?:', np.sum(y[i]))

    #CDF
    partial_sum = 0
    for j in range (0, b):
        print(b-j)
        partial_sum = partial_sum + y[i][j]
        cdf_y[i][j] = partial_sum

    # GRAPHIC
    # ax1.plot(np.arange(1,len(y[i][:b])+1),y[i][:b], linewidth = 0.5)#, label='PMF_'+str(c[i])+' '+ Adata[1][c[i]] +', '+str(b)+' words.')
    ax1.plot(np.arange(1, len(y[i][:b]) + 1), cdf_y[i][:b],
             linewidth=0.5, label='CDF_' + str(c[i]) + ' ' + Adata[1][c[i]] + ', ' + str(b) + ' words.')

grid = np.arange(S, n_counties)
couples = list(itertools.permutations(np.arange(S, n_counties), 2))
D = np.zeros((len(grid), len(grid)))
p = np.zeros((len(grid), len(grid)))
myD = np.zeros((len(grid), len(grid)))

'''
# Filling triangular matrix
for i in range(0, len(grid)):
    for j in range(i + 1, len(grid)):
        if cdf_y[i + S] > cdf_y[j + S] :
            abs( )
        else:
            ascii()
    print(len(grid) - i)
combinations = len(grid) * (len(grid) - 1) * 0.5
'''

#lent =[len(cdf_y[l])for l in range (n_counties) ]
#print(lent)


# Setting log-log scale
ax1.set_yscale('log')
ax1.set_xscale('log')

# Axis titles
ax1.set_xlabel('Ranking')
ax1.set_ylabel('Normalized frequency')

ax1.legend(loc='best')

#print(len(cdf_y[1+S][:int(min(word[1+S],word[4+S]))]))
#print(len(cdf_y[4+S][:int(min(word[1+S],word[4+S]))]))

#Filling triangular matrix
for i in range (0, len(grid)):
    for j in range (i+1, len(grid)):
        D[i][j], p[i][j] = stats.ks_2samp(cdf_y[i+S],cdf_y[j+S])
        myD[i][j] = abs(np.amax(   cdf_y[i+S][:int( min( word[i+S],word[j+S] ) )] - cdf_y[j+S][:int ( min( word[i+S],word[j+S] )) ]))

figDD = plt.figure()
axD = figDD.add_subplot(1, 1, 1)
#print(, len(g_D))
axD.scatter(D,myD, s = 0.5)
# Axis titles
axD.set_xlabel('KS ')
axD.set_ylabel('my KS')
print(cdf_y[1])
print(cdf_y[2])

plt.show()
