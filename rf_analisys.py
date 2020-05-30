import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)

#total no. of counties
N_counties = len(Adata[1][:]) -1

#choose the county by a number between [1,N_counties]
c = np.random.choice(np.arange(1, N_counties + 1 ))
c = 746
county = Adata[1][c]

#chooose the number of words you need; default is the maximum
word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False, index_col=[0])
word = int(word.iloc[c, 1])
print(word)

#IMPORTING DATA
filename = r"C:\Users\Salvo\Desktop\IFISC\data\r_f_county_w\freq_ranking_"+county+'(' + str(c) +')_'+str(word)+"w.csv"
dat = pd.read_csv(filename)
print(dat)

#ZIPF SAMPLE, alpha
alpha = 1.1877957588506238
#alpha = 1 + 1/alpha
sizer = word*1000


#first sample

#random_zipf= np.random.zipf(alpha, size=sizer)
#random_zipf = np.int64(- np.sort(-random_zipf))
#norm_fact = np.sum(random_zipf)
#print(norm_fact)
#random_zipf = random_zipf/norm_fact
#print('Are the 1# sample data normalized?',np.sum(random_zipf))
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

#second sample
random_zipf1 = stats.zipf.sf(np.arange(1,word+1),a = alpha)
random_zipf2 = stats.zipf.pmf(np.arange(1,word+1),a = alpha)
x_sample = np.arange(1,len(random_zipf1)+1)
#random_zipf1 = np.int64 (- np.sort(-random_zipf1))
#norm_fact1 = np.sum(random_zipf1)
#print(norm_fact1)
#random_zipf1 = random_zipf1/np.sum(random_zipf1)
print('Are the 2# sample data normalized?',np.sum(random_zipf1))


#Real normalized data

#PDF
x = np.arange(1, word + 1)
y = np.int64(dat.iloc[:, 2])
sum_y = np.sum(y)
y = y/sum_y
norm_sum_y = np.sum(y)
print('Are the data normalized ?:', norm_sum_y)

#Survival function
sf_y = np.zeros(word)
for i in range (0, word-1):
    sf_y[i] = norm_sum_y - np.sum(y[i+1:word])
    print(sf_y[i])
sf_y[word-1] =  norm_sum_y
sf_y = 1 - sf_y
#sf_y = sf_y/np.sum(sf_y)

#parole = dat.iloc[:, 1]
#print(x, y, parole)


# GRAPHIC
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)


#x_mysample = np.arange(1,int(len(random_zipf))+1)
#print(len(x_mysample))
ax1.scatter(x, y, label=county, color='red', s = 1)
ax1.scatter(x, sf_y, label='survival function: '+county, color='xkcd:neon purple', s = 1)
#ax1.scatter(x_mysample,random_zipf, label="Zipf's sampled data", s=1, color = 'grey')
ax1.scatter(x_sample,random_zipf1,  label="Zipf's scipy sampled data, SF", s=1, color='orange')
ax1.scatter(x_sample,random_zipf2,  label="Zipf's scipy sampled data, PMF", s=1, color='black')

#ax1.grid(color='xkcd:cherry', linestyle='--', which='major', zorder=1, alpha=0.75)
#ax1.set_xticks(np.arange(0, word + 1), minor=False)
#for i, parole in enumerate(parole):
    #ax1.annotate(parole, (x[i] - 0.2, y[i] + 1.2))



#LINEAR REGRESSION

#a = starting point, b = end point
a = 20
b = 1000

#LR_PDF on real data
x = np.arange(a, b)
y = y[a:b]
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x), np.log10(y))
print('LR_PDF m for real data =', slope)
print('LR_PDF q  for real data =', intercept)

#LR_SF on real data
x = np.arange(a, b)
sf_y = sf_y[a:b]
slope_sf, intercept_sf, r_value_sf, p_value_sf, std_err_sf = stats.linregress(np.log10(x), np.log10(sf_y))
print('LR_SF m for real data =', slope_sf)
print('LR_SF q for real data =', intercept_sf)


#plotting LR_PMF for real data
c1 = 10**intercept
x1 = np.arange(1, word+1)
y1 = (x1**slope) * c1
#print(np.sum(y1))
ax1.plot(x1, y1, label="LR_PDF on real data", linewidth=1, color='green')

#plotting LR_SF for real data
c1_sf = 10**intercept_sf
x1_sf = np.arange(1, word+1)
y1_sf = (x1_sf**slope_sf) * c1_sf
#print(np.sum(y1))
ax1.plot(x1_sf, y1_sf, label="LR_SF on real data", linewidth=1, color='xkcd:yellow brown')



#LR on sampled data
#slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(np.log10(x), np.log10(random_zipf[a:b] ))
#print('LR slope for sampled data =', slope1)
#print('LR intercept for sampled data =', intercept1)

slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(np.log10(x), np.log10(random_zipf1[a:b] ))
print('LR_SF m for sampled data = ', slope2)
print('LR_SF q for sampled data = ', intercept2)
slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(np.log10(x), np.log10(random_zipf1[a:b] ))
print('LR_PMF m for sampled data = ', slope3)
print('LR_PMF q for sampled data = ', intercept3)



#Plotting sampled data and LR on sampled data, first sample
#y2 = x1**slope1 * 10**intercept1
#ax1.plot(x1,y2 , linewidth=1, color='purple', label = "LR on Zipf's sampled data")

#Plotting sampled data and LR on sampled data, second sample SF
y3 = x1**slope2 * 10**intercept2
ax1.plot(x1,y3 , linewidth=1, color='blue', label = "LR on Zipf's sampled data 2, SF")

#Plotting sampled data and LR on sampled data, second sample PMF
y4 = x1**slope3 * 10**intercept3
ax1.plot(x1,y4 , linewidth=1, color='brown', label = "LR on Zipf's sampled data 2, PMF")


#Setting log-log scale
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend(loc='best')


ax1.set_title('Real and sampled scattered data with their linear regression plotted')

ax1.set_xlabel('Rank')
ax1.set_ylabel('Normalized frequency')

plt.savefig('C:\\Users\Salvo\Desktop\IFISC\data\graphics\LR_test_sampled&real_'+str(county)+'_('+str(c)+').svg')

plt.show()










#x1 = np.arange(a, b, 0.01)
#y_sample = (1/2.612) * (x1**-1.12)
#y_zipf_sample =
#print(y_sample[0],x[0])
#ax1.plot(x1,y_sample, label="l-r on Zipf's sample",color='orange', linewidth=1)

'''
sloppete = np.zeros(b-a+1)
interceptete = np.zeros(b-a+1)
std_errete = np.zeros(b-a+1)

sloppete1 = np.zeros(b-a+1)
interceptete1 = np.zeros(b-a+1)
std_errete1 = np.zeros(b-a+1)
'''

'''
sloppete = np.zeros(b-a)
interceptete = np.zeros(b-a)
std_errete = np.zeros(b-a)

for i in range (a, b, 1) :
    x = np.arange(i, word)

    print(x,y)
    print(b -i)
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x),np.log10(y.drop(y.index[np.arange(0,i)], axis = 0)))

    sloppete[i-a] = slope
    interceptete[i-a] = intercept
    std_errete[i-a] =std_err
'''



'''
for i in range(a, b, 1):
    x = np.arange(a, b+a-i)

    print(x, y)
    print(b - i)
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x), np.log10(y.iloc[a:(b+a-i)]))
    slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(np.log10(x), np.log10(zipf_sample[a:(b+a-i)]))

    sloppete[b-i] = slope
    interceptete[b-i] = intercept
    std_errete[b-i] = std_err

    sloppete1[b - i] = slope1
    interceptete1[b - i] = intercept1
    std_errete1[b - i] = std_err1


    #x1 = np.arange(1, word + 1, 0.1)
    #y1 = 10**intercept*x1**(slope)
    #ax1.plot(x1,y1)

    #print(slope,10**(intercept), std_err)
    #print(parole)

x = np.arange(a, b+1)
#x = np.arange(a, b, 1)
ax2 = fig.add_subplot(1, 1, 1)
ax2.scatter(x, sloppete, s = 0.5 )
ax2.scatter(x_sample, sloppete1, s = 0.5 )
#ax2.scatter(x, std_errete,s = 0.5)
#filename = r"C:\\Users\Salvo\Desktop\IFISC\data\graphics\freq_ranking_"+county+'(' + str(c) +')_'+str(word)+"w.eps"
#fig.savefig('C:\\Users\Salvo\Desktop\IFISC\data\graphics', format='eps')
plt.show()

'''


'''
#Generating random Zipf distributed data
def randzipf(x_min, x_max, a):
    x = np.random.rand(x_max-x_min)
    C = (a+1) / (x_max**(a+1)- x_min**(a+1))
    y = C/(a+1) * (x**(a+1)- x_min**(a+1))
    X =  ((x_max**(a+1)- x_min**(a+1))*y + x_min**(a+1))**(1/(a+1))
    X = -np.sort(-X)
    X = X/(np.sum(X))
    return X
'''



