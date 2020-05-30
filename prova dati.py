import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from plfit import plfit
from itertools import product
import pickle
import csv
from sklearn.externals import joblib


a = pd.read_csv(r'C:\Users\Salvo\Desktop\da inviare\js2D_value.csv', sep=",", low_memory=False)
print(a)
print(a.iloc[2,3], a.iloc[3,2])
print(a.iloc[2, :])
print(a.iloc[:,2])


ks_file = np.load(r'C:\Users\Salvo\Desktop\da inviare\ks2D_value_4723201_combinations.npz')
js_file = np.load(r'C:\Users\Salvo\Desktop\da inviare\js2D_value_4723201_combinations.npz')

print(ks_file.files)
print(js_file.files)


ks_D = ks_file['ks_D']
js_D = js_file['js_D']

print(ks_D)
print(js_D)

with open(r'C:\Users\Salvo\Desktop\da inviare\ks2D_value.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(
        ks_D
    )
print('CAZZI IN CULO')
with open(r'C:\Users\Salvo\Desktop\da inviare\js2D_value.csv', 'w', newline='') as csvfile1:
    writer = csv.writer(csvfile1)
    writer.writerows(js_D)
print('MINCHIA TESA')


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
'''
#JS
combinations = 4723201
#Importing JS
js_distance = np.load(r'C:\\Users\Salvo\Desktop\IFISC\data\all_counties\js\run\jsD_value_'+str(int(combinations))+'_combinations.npz')
print(js_distance.files)

#Extracting arrays from file
js_D = js_distance['jd_D']
js_grid = js_distance['counties']







'''
# List of total number of words per county
tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,index_col=[0])
print(tot_word)



#with open(r'C:\Users\Salvo\Desktop\y2.pkl', 'rb') as f:
 #   y = pickle.load(f)
#print(y[0], y[3074])

'''
# IMPORTING DATA
filename = "C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\\row_frequencies\\frequencies_connecticut;windham_(3075).txt"
dat = pd.read_csv(filename, sep=" ", header=None)
ultima = np.asarray(dat.iloc[:, 0].to_numpy(), dtype=float)
print(len(ultima),ultima)
ultima = dat.iloc[:, 0]
ultima = - np.sort(-ultima)
ultima = ultima / np.sum(ultima)
y.append(ultima)
'''


#print(y[3074], len(y[3074]))

print('ciao ora inizia il for I')
cdf_y = [[]]
# CDF
for i in range(1, 3075):
    b = int(tot_word.iloc[i, 1])


print(cdf_y[0])#, cdf_y[1999])

print('Ho calcolato le cdf')


#np.savez_compressed(r'C:\Users\Salvo\Desktop\IFISC\data\cdf1_y', cdf_y = cdf_y)



#with open(r'C:\Users\Salvo\Desktop\IFISC\data\y2.pkl', 'wb') as f:
 #   pickle.dump(y, f)
#filename = r'C:\Users\Salvo\Desktop\IFISC\data\y2.csv'









#KS
#Importing KS
combinations = 4723201
myks_file = np.load(r'C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\ks\\run\myD_value_'+str(int(combinations))+'_combinations.npz')
#Extracting arrays from file
myD = myks_file['D']
ks_grid = myks_file['counties']


values = np.zeros((3074, 1))
new_row = np.zeros((1, 3075))
cdf_y = np.zeros(70000)
cdf_x = np.zeros(39898)
partial_sum = 0

for j in range(39898):
    partial_sum = partial_sum + y[3074][j]
    cdf_x[j] = partial_sum

print('ecco la puttana', cdf_x)

for i in range (1,3075):
    b = int(tot_word.iloc[i, 1])

    partial_sum = 0
    # print(b, len(y[i]))

    for j in range(b):
        partial_sum = partial_sum + y[i - 1][j]
        cdf_y[j] = partial_sum
    print(3075 - i)


    minNWords = int(min(b, 39898))
    values[i-1] = np.amax(abs(np.array(cdf_y[:minNWords]) - np.array(cdf_x[:minNWords])))
#values[3074]= 0
print(values)
print(myD.shape)

myD = np.hstack((myD,values))
myD = np.vstack([myD, new_row])

print(myD, myD.shape)
#np.savez(r'C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\js\\run\ks2D_value_'+str(int(combinations))+'_combinations', ks_D=myD, counties = ks_grid)

#, js_D[0][3074], js_D[3074][3074])
#dat = pd.read_csv(filename, sep=" ", header=None)
#print(dat)

#with open(r'C:\Users\Salvo\Desktop\IFISC\data\y2.csv', 'w', newline='') as csvfile:
#    writer = csv.writer(csvfile)
#    writer.writerows(y)




'''

with open(r'C:\\Users\Salvo\Desktop\IFISC\data\cdf1_y.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(cdf_y)
    print('fattoh')



with open(r'C:\\Users\Salvo\Desktop\IFISC\data\cdf1_y.pkl', 'wb') as f:
    pickle.dump(cdf_y, f)

#KS
combinations = 4723201
#myks_file = np.load(r'C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\ks\\run\myD_value_'+str(int(combinations))+'_combinations.npz')
myD = myks_file['D']
ks_grid = myks_file['counties']
print(myD, len(myD[0][:]))


N_ks = np.zeros((len(ks_grid),(len(ks_grid) )))
#js_Dfull = np.ones((len(js_grid),(len(js_grid) )))
ks_Dfull = myD
for i in range (0, len(ks_grid)):
    for j in range (i+1, len(ks_grid)):
        ks_Dfull[j][i] = myD[i][j]

print(myD)
print('AAAAAAAAAAAAAAAAAAAAAA')
print(ks_Dfull)


for k in range (10):
    N_grid = np.arange(3074)
    np.random.shuffle(N_grid)
    print(N_grid)

    for i in range(len(ks_grid)):
        for j in range (i+1,len(ks_grid)):
            #print(N_grid[i])
            #print(N_grid[j])
            #print(js_Dfull[N_grid[i]][N_grid[j]])

            N_ks[i][j] = ks_Dfull[N_grid[i]][N_grid[j]]
    print(N_ks)
    np.savez(r'C:\\Users\Salvo\Desktop\IFISC\data\\New_Dist\\N_ks_D_'+str(k), N_ks_D = N_ks, N_grid = N_grid)


#JS
combinations = 4723201
js_distance = np.load(r'C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\js\\run\jsD_value_'+str(int(combinations))+'_combinations.npz')
print(js_distance.files)
js_D = js_distance['jd_D']
js_grid = js_distance['counties']

N_js = np.zeros((len(js_grid),(len(js_grid) )))
#js_Dfull = np.ones((len(js_grid),(len(js_grid) )))
js_Dfull = js_D
for i in range (0, len(js_grid)):
    for j in range (i+1, len(js_grid)):
        js_Dfull[j][i] = js_D[i][j]
print(js_D)
print('AAAAAAAAAAAAAAAAAAAAAA')
print(js_Dfull)


for k in range (10):
    N_grid = np.arange(3074)
    np.random.shuffle(N_grid)
    print(N_grid)

    for i in range(len(js_grid)):
        for j in range (i+1,len(js_grid)):
            #print(N_grid[i])
            #print(N_grid[j])
            #print(js_Dfull[N_grid[i]][N_grid[j]])

            N_js[i][j] = js_Dfull[N_grid[i]][N_grid[j]]
    print(N_js)
    np.savez(r'C:\\Users\Salvo\Desktop\IFISC\data\\New_Dist\\N_js_D_'+str(k), N_js_D = N_js, N_grid = N_grid)


print('fatto boxhy')


# List of total number of words per county
tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,index_col=[0])
#print(tot_word)

'''

'''
Adata = pd.read_csv(r'C:\\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)
print(Adata)

county_list = Adata.iloc[:,1].to_list()
lista1 = []
lista2 = []
print(county_list)

indexes = []

for i in range (1,len(county_list)):
    indexes.extend([pos for pos, char in enumerate(county_list[i]) if char == ';'])
    print(indexes[i-1])
    lista1.append(county_list[i][ : indexes[i-1]])
    lista2.append(county_list[i][  indexes[i-1]+1:])



print(lista1)

lista1 = sorted(set(lista1))
lista2 = sorted(set(lista2))


#total no. of counties
N_counties = len(Adata[1][:]) -1

#choose the county by a number between [1,N_counties]
c = np.random.choice(np.arange(1, N_counties + 1 ))
c = 3000
county = Adata[1][c]

#chooose the number of words you need; default is the maximum
word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False, index_col=[0])
word = int(word.iloc[c, 1])
print(word)

#IMPORTING DATA
filename = r"C:\\Users\Salvo\Desktop\IFISC\data\\r_f_county_w\\freq_ranking_"+county+'(' + str(c) +')_'+str(word)+"w.csv"
dat = pd.read_csv(filename)
print(dat)

a = np.arange(10)
b = np.append(a, np.zeros(6))
print('b: ', b)

'''

'''
#PDF
x = np.arange(1, word + 1)
y = np.int64(dat.iloc[:, 2])
sum_y = np.sum(y)
#y = y/sum_y
#norm_sum_y = np.sum(y)
#print('Are the data normalized ?:', norm_sum_y)
for i in range (0, len(y)):
    print('y[i]',y[i])

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



word = 5000
a = 3000
x_min = 1
x = np.arange(x_min, word)

Alpha = np.zeros(a)
alpha = np.zeros(a)


for x_min in range (1, a):
    Alpha[x_min] = discrete_alpha_mle(x, xmin = x_min)



sum_ln_x_i = 0
n = word - x_min
x = np.arange(x_min, word)
for i in range(1, n):
    sum_ln_x_i = sum_ln_x_i + np.log(x[i])

for x_min in range (1, a):
    n = word - x_min
    print(a-x_min)
    sum_ln_x_i = sum_ln_x_i - np.log(x_min)
    #print(np.log(1/(x_min - 0.5)), sum_ln_x_i)
    alpha[x_min-1] = 1 + n*(sum_ln_x_i - n*np.log(1/(x_min))  )**-1
    #print(alpha)
print(len(alpha))

fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

ax1.scatter(np.arange(1,a+1), alpha,label = 'my estimator', s = 1)
ax1.scatter(np.arange(1,a+1), Alpha,label = 'package estimator', s = 1)

ax1.set_xlabel('x_min')
ax1.set_ylabel('alpha')
ax1.set_title('MLE alpha estimator, alpha vs x_min (using ranking)')
ax1.set_xscale('log')
ax1.set_yscale('linear')


sigma = (alpha - 1)/np.sqrt(n)
print('sigma',sigma)
print('alpha',Alpha)

X = np.arange(1, 8000)
MyPL = plfit(X)

plt.show()
'''
'''
X = np.arange(50, 330) #x_min
Y = np.arange(x_min, n, 25) #n
X, Y = np.meshgrid(X, Y)
alpha = 1 + Y*(sum_ln_x_i - np.log(1/(X - 0.5)))**-1
alpha_sample = 1 + Y*(sum_ln_x_i - np.log(1/(X - 0.5)))**-1
sigma = (alpha - 1)/np.sqrt(n)
print(alpha)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, alpha)
#ax.plot_surface(X, Y, sigma)

plt.show()
'''
'''
alpha = np.zeros(1)

for i in range (1,n):
    alpha = alpha + (np.log(1/(x_min - 0.5)))

alpha = 1 + n*(sum_ln_x_i -alpha)**-1
sigma = (alpha - 1)/np.sqrt(n)
print(alpha, sigma)


alpha = np.zeros(n-1000)
for k in range (1, n-1000):
    x_min = k
    print(k)
    for i in range (1,n):
        alpha[k] = alpha[k] + (np.log(1/(x_min - 0.5)))

    alpha[k] = 1 + n*(sum_ln_x_i -alpha[k])**-1
sigma = (alpha - 1)/np.sqrt(n)
print(alpha, sigma)

x = np.arange(1, n-999)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.scatter (x, alpha, s = 0.5)

plt.show()

'''



'''
for k in range(50, 330) :
    x_min = k
    print('k:', 300 - k)
    for j in range (1000, 8000, 25):
        n = j
        a = 0
        print('j:', 8000-j)
        for i in range(1,n-x_min):
            x = np.arange(x_min, n + 1)
            a = a + (np.log(x[i]/(x_min - 0.5)))
            Z[] =

        alpha = 1 + n*a**-1
        print(alpha)
'''







'''
Adata = pd.read_csv(r'C:\\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)
prova = pd.DataFrame()



which_county = 1000
N_counties = Adata[1][:]
county = Adata[1][1]

print(N_counties, len(N_counties))


word = 20
c = 3
z = np.zeros(word)
z[0] = int(c)
rfcw = pd.DataFrame({'z': z})

print(z, rfcw)




#Bdata = pd.read_csv(r'C:\\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_B_PBW.txt', sep=",", header= None, low_memory=False)
#Cdata = pd.read_csv(r'C:\\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_C_PBW.txt', sep=",", header= None, low_memory=False)




#ATdata = np.transpose(Adata)
#BTdata = np.transpose(Bdata)
#CTdata = np.transpose(Cdata)
#print(ATdata)

#for i in range (1,3801):
    #print (ATdata[0][i])
#print(ATdata)

#N_counties = ATdata.shape[1] -1
#words = ATdata.shape[0] - 4
#ount = words
#print(N_counties)

#counties = np.zeros((3,words,count))


#for i in range (1, 4):
 #   for k in range (4,words):
  #      counties[i][k][k] = [[ATdata[1][i]][ATdata[k][0]][ATdata[k][i]]]
#print(counties)
'''
'''
Adata = Adata.drop(columns=[0, 2, 3])
print(Adata)
Adata = Adata.iloc[:2, :]
Adata = Adata.T
Adata = Adata.drop(index=[1])
print(Adata)
Adata[1] = pd.to_numeric(Adata[1])
Adata = Adata.sort_values(by = [1], ascending=False)
Adata = Adata.reset_index()

print(Adata)
'''

'''
random_zipf = np.asarray(random_zipf, dtype=np.float)
print('MAAAAAX', np.max(random_zipf))
binss = np.arange(1,np.max(random_zipf))
random_zipf = np.histogram(random_zipf, bins = binss)
random_zipf = np.array(random_zipf)
random_zipf_a = np.array(random_zipf[0])
random_zipf_b = np.array(random_zipf[1])
#print(y, y0, y0[35], y1)
no0_idx = np.array(np.nonzero(random_zipf_b))
print('AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', no0_idx[0])
no0_idx = no0_idx[0]
random_zipf_b = -np.sort(-random_zipf_b)
random_zipf = random_zipf_b[no0_idx]
#random_zipf = random_zipf[np.nonzero(random_zipf)[0]]
random_zipf = random_zipf/np.sum(random_zipf)
'''