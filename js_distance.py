import numpy as np
import pandas as pd
import pickle

#kl divergence
def kl_D (p,q):

    kl_d = np.sum(p * np.log(p / q))
    #print('1:',kl_d)
    if np.isnan(kl_d)  :
        kl_d = 0
    #print('2:',kl_d)
    return kl_d

#js metric, This needs the kl divergence function
def js (p, q):
    #checking leght
    a = len(p)-len(q)
    #This if elif appends len(p)-len(q) zeros to the shortest distribution
    #So they are the same lenght
    if a <0 :
        p = np.append(p, np.zeros(-a))
    elif a>0:
        q = np.append(q, np.zeros(a))

    #m = mean p,q
    m = (p+q)*0.5

    #js distance
    js_d = np.sqrt(0.5* (kl_D(p,m)+kl_D(q,m)) )

    return js_d

Adata = pd.read_csv(r'C:\\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header=None,
                    low_memory=False)

# total no. of counties
TOT_counties = len(Adata[1][:]) - 1
print(TOT_counties)

# List of total number of words per county
tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,
                       index_col=[0])
print(tot_word)

which_county = np.arange(1,TOT_counties)
# number of words of the i-th county
word = np.zeros(len(which_county))
# indx o i-th county
c = which_county
print('ciao')


which_county = np.arange(1,TOT_counties-1 +1)

with open(r'C:\Users\Salvo\Desktop\y.pkl', 'rb') as f:
    y = pickle.load(f)

grid = which_county
jd_D = np.zeros((len(which_county), len(which_county)))
combinations = len(grid)*(len(grid)-1)*0.5

"""
Probably the filling doesn't include the last county.
I added it manually afterward. This should be fixed but the distance calulation for
a single couple of counties is ok
"""
#Filling upper triangular matrix
for i in range (0, len(grid)):
    #print(i)
    for j in range (i+1, len(grid)):
        #print(jd_D[i][j])
        jd_D[i][j]= js(y[i ], y[j ])
        #np.savez(r"C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\js\jsD_value_(" + str(i) + '-' + str(j) + ")", counties=str(i)+','+str(j), jd_D=jd_D[i][j])

    print(len(grid)-i)
combinations = len(grid)*(len(grid)-1)*0.5
np.savez(r'C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\js\\run\jsD_value_'+str(int(combinations))+'_combinations', jd_D=jd_D, counties = grid)

