import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import scipy.stats as stats



Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)
print(Adata)
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
y = []
cdf_y = []

#Initializing variables
x = np.arange(1, 70000)
for i in range (S,n_counties):
    word[i] = int(tot_word.iloc[c[i], 1])
    b = int(word[i])
    y.append([])
    cdf_y.append([])
print(cdf_y)
#print(y)

#GRAPHIC
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)


for i in range (S,n_counties):
    word[i] = int(tot_word.iloc[c[i], 1])
    #print(word[i])
    #print(Adata[1][c[i]])
    #print(y[i])

    #IMPORTING DATA
    filename = "C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\\row_frequencies\\frequencies_" + str(Adata[1][i]) + "_(" + str(c[i]) + ")" + ".txt"
    dat = pd.read_csv(filename, sep=" ", header=None)
    #print(dat)

    b = int(word[i])

    #Real normalized data
    #PMF
    #y[i] = dat.iloc[:,0]
    #y[i] = - np.sort(-y[i])
    #y[i] = y[i]/np.sum(y[i])
    #print('Are the data normalized ?:', np.sum(y[i]))

    #CDF
    partial_sum = 0
    '''
    for j in range (0, b):

        partial_sum = partial_sum + y[i][j]
        cdf_y[i].append(partial_sum)

    print(cdf_y[i])
    f = open("C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\\cdf\\cdf_" + str(
        Adata[1][i]) + "_(" + str(c[i]) + ")" + ".txt", "a+")
    for j in range(len(cdf_y[i])):
        f.write(str(cdf_y[i][j]) + '\n')
    f.close()
    print(n_counties-i)
'''
print('WORDI I:', word)

np.savez(r"C:\Users\Salvo\Desktop\IFISC\data\all_counties\cdf\cdf_(" + str(1) + '-' + str(TOT_counties) + ")", cdf_y=cdf_y)


