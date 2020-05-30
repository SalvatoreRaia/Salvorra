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
c = np.arange(0, n_counties)
print('ciao')

#Initializing variables
x = np.arange(1, 70000)
for i in range (S,n_counties):
    word[i] = int(tot_word.iloc[c[i], 1])
    b = int(word[i])

    y = [[0 for indiceacaso in range(1)] for indiceacaso2 in range(n_counties)]
    #cdf_y = [[0 for indiceacaso in range(1)] for indiceacaso2 in range(n_counties)]

#print(y)

#GRAPHIC
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

cdf_y = []
for i in range (S,n_counties):
    word[i-S] = int(tot_word.iloc[c[i], 1])
    print(word[i])
    print(Adata[1][c[i]])
    print(y[i])

    #IMPORTING DATA
    filename = "C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\\row_frequencies\\frequencies_" + str(Adata[1][i]) + "_(" + str(c[i]) + ")" + ".txt"
    dat = pd.read_csv(filename, sep=" ", header=None)
    #print(dat)
    row = []
    b = int(word[i-S])

    #Real normalized data
    #PMF
    y[i] = dat.iloc[:,0]
    y[i] = - np.sort(-y[i])
    y[i] = y[i]/np.sum(y[i])
    print('Are the data normalized ?:', np.sum(y[i]))

    #CDF
    #print('primo valore y:', y[i][b-1])
    partial_sum = [0.]
    for j in range (0, b):
        #print(partial_sum[0])
        partial_sum[0] = partial_sum[0] + y[i][j]
        row.append(float(partial_sum[0]))
        #print(row)
    #print(row)

    cdf_y.append(row)
    print(cdf_y)
    print('SUcchia cazzi',cdf_y[i-S])

    #GRAPHIC
    #ax1.plot(np.arange(1,len(y[i][:b])+1),y[i][:b], linewidth = 0.5)#, label='PMF_'+str(c[i])+' '+ Adata[1][c[i]] +', '+str(b)+' words.')
    #ax1.plot(np.arange(1, len(y[i][:b])+1 ), cdf_y[i][:b], linewidth=0.5)#,label='CDF_'+str(c[i]) + ' ' + Adata[1][c[i]] + ', ' + str(b) + ' words.')
print('la prima',cdf_y[0])
print('lultimo valore',cdf_y[3][-1:])
#KS Statistics
#compressed_cdf_y =
grid = np.arange(S, n_counties)

couples = list(itertools.permutations(np.arange(S, n_counties), 2))
D = np.zeros((len(grid), len(grid)))
myD = np.zeros((len(grid), len(grid)))
otherD = np.zeros((len(grid), len(grid)))

p = np.zeros((len(grid), len(grid)))
#print(D)
ok_couples = []
ok_couples_pi = []

print('la prima y', y[0])
print('Le coppie',couples)
print('check lunghezze :', len(cdf_y[0]), word[0])

#Filling triangular matrix
for i in range (0, len(grid)):
    for j in range (i+1, len(grid)):
        #lengti = int(min(len(cdf_y[i]), len(cdf_y[j])))
        #unae = cdf_y[i+S][:lengti]
        #ddue = cdf_y[j+S][:lengti]
        #D[i][j], p[i][j] = stats.ks_2samp(cdf_y[i+S],cdf_y[j+S])
        #otherD[i][j] = abs(np.amax(np.array(unae) - np.array(ddue)))
        myD[i][j] = np.amax(abs(np.array(cdf_y[i][:int(min(word[i], word[j]))]) - np.array(cdf_y[j][:int(min(word[i], word[j]))])))
        #print(i, word[j+S], word[i+S])
        np.savez(r"C:\Users\Salvo\Desktop\IFISC\data\all_counties\ks\myD_value_(" + str(i) + '-' + str(j) + ")", counties=str(i)+','+str(j), myD=myD[i][j])

    print(len(grid)-i)
combinations = len(grid)*(len(grid)-1)*0.5
np.savez(r'C:\Users\Salvo\Desktop\IFISC\data\all_counties\ks\run\myD_value_t'+str(int(combinations))+'_combinations', D=myD, counties = grid)


print('WORDI I:', word)
#print(ok_couples)
#print(ok_couples_pi)
#print(len(ok_couples))

#differenza = myD - otherD
#print('ORA PRInto')
#for i in range(len(differenza)):
 #   print(differenza[i])

#ax1.scatter(np.arange())


#Figure settings

#Setting log-log scale
ax1.set_yscale('log')
ax1.set_xscale('log')

#Axis titles
ax1.set_xlabel('Ranking')
ax1.set_ylabel('Normalized frequency')

#Axis grid
#ax1.grid(color='xkcd:cherry', linestyle='--', which='major', zorder=1, alpha=0.75)
#ax1.set_xticks(np.arange(0, 70000 + 1), minor=False)

ax1.legend(loc='best')
#ax1.scatter(myD,otherD, s = 0.5, color = 'red')



#plt.show()
#ax1.plot(np.arange(len(sf_y)),sf_y, linewidth = 0.5, label='CDF')






