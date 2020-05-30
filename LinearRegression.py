import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from numpy import sum


######################################
"""
This program plot a single distribution
1) calculate double linear regression (LR) on it
2) Use MLE to calculate ZIpf alpha parameter and with prints on screen, you can compare it with the LR alpha
3) Comparison with sampled data using scipy func. On sampled data it does also linear regression to verify that the alpha
resulting from LR on sampled data , is the same of the one used to generate the data


You can rewrite the code to support this analysis for more counties at the same time just importing the file in which all
the counties are stored and looping this program for the counties you want to analyse
"""
######################################


def sigma(alpha, n):
    """
    Clauset et al 2007 equation 3.2:
        sigma = (alpha-1)/sqrt(n)
    """
    return (alpha-1.) / n**0.5


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


def discrete_ksD(data, xmin, alpha):
    """
    given a sorted data set, a minimum, and an alpha, returns the power law ks-test
    D value w/data

    The returned value is the "D" parameter in the ks test

    (this is implemented differently from the continuous version because there
    are potentially multiple identical points that need comparison to the power
    law)
    """

    zz = np.sort(data[data >= xmin])
    nn = float(len(zz))
    if nn < 2:
        return np.inf
    # cx = np.arange(nn,dtype='float')/float(nn)
    # cf = 1.0-(zz/xmin)**(1.0-alpha)
    model_cdf = 1.0 - (zz.astype('float') / float(xmin)) ** (1.0 - alpha)
    data_cdf = np.searchsorted(zz, zz, side='left') / (float(nn))
    ks = max(abs(data_cdf - model_cdf))
    return ks



Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)

#total no. of counties
N_counties = len(Adata[1][:]) -1

#choose the county by a number between [1,N_counties]
#c = np.random.choice(np.arange(1, N_counties + 1 ))
#print(c)
c = 1968
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
alpha = 1.8976294710967334
#alpha = 1 + 1/alpha
sizer = word*10000


#first sample

#random_zipf= np.random.zipf(alpha, size=sizer)
#random_zipf = np.int64(- np.sort(-random_zipf))
#norm_fact = np.sum(random_zipf)
#print(norm_fact)
#random_zipf = random_zipf/norm_fact
#print('Are the 1# sample data normalized?',np.sum(random_zipf))


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

#PMF (Probability mass function, eqeuivalent for PDF in descrete cases)
x = np.arange(1, word + 1)
y = np.int64(dat.iloc[:, 2])
sum_y = np.sum(y)
#Normalizing
y = y/sum_y
norm_sum_y = np.sum(y)
print('Are the data normalized ?:', norm_sum_y)


# GRAPHIC
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

#x_mysample = np.arange(1,int(len(random_zipf))+1)
#print(len(x_mysample))
ax1.scatter(x, y, label=county, color='purple', s = 0.5)

#Plotting sampled data
#ax1.scatter(x, sf_y, label='survival function: '+county, color='xkcd:neon purple', s = 0.1)
#ax1.scatter(x_mysample,random_zipf, label="Zipf's sampled data", s=1, color = 'grey')
#ax1.scatter(x_sample,random_zipf1,  label="Zipf's scipy sampled data, SF", s=0.1, color='black')
#ax1.scatter(x_sample,random_zipf2,  label="Zipf's scipy sampled data, PMF", s=2.5, color='pink')



#axis settings
#ax1.grid(color='xkcd:cherry', linestyle='--', which='major', zorder=1, alpha=0.75)
#ax1.set_xticks(np.arange(0, word + 1), minor=False)




#LINEAR REGRESSION (LR)

#a = starting point, b = end point
a = 1
b = 2
#LR_PDF on real data
x = np.arange(a, b)
ys = y[a:b]
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x), np.log10(ys))
print('LR_PDF m for real data =', slope)
print('LR_PDF q  for real data =', intercept)
b1 = word

#plotting dot of behaviour separation
ax1.scatter(x=b, y=(b**slope) * 10**intercept, s = 15, color = 'black', label = 'behave separetor', zorder=500)

xr = np.arange(b, b1)
yr = y[b:b1]
sloper, interceptr, r_valuer, p_valuer, std_errr = stats.linregress(np.log10(xr), np.log10(yr))
print('LR_PDF m for real data =', sloper)
print('LR_PDF q  for real data =', interceptr)


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

print('MLE_alpha, behaviour 1 :',discrete_alpha_mle(np.arange(1, word), a))
print('MLE_alpha, behaviour 2 :',discrete_alpha_mle(np.arange(1, word), b))
#print('MLE_alpha, behaviour_all_words, xmin graphically :',discrete_alpha_mle(np.arange(a, word), a))


#plotting LR_PMF for real data


#1 behaviour
c1 = 10**intercept
x1 = np.arange(1, b+1)
y1 = (x1**slope) * c1
ax1.plot(x1, y1, label="LR_PMF on real data, $\\alpha_1$ = " +str(slope), linewidth=1, color='red', zorder=100)

#2 behaviour
cr = 10**interceptr
xr = np.arange(b, word)
yr = (xr**sloper) * cr
#print(np.sum(y1))
ax1.plot(xr, yr, label="LR_PMF on real data, $\\alpha_2$ = "+str(sloper), linewidth=1, color='orange', zorder=120)





#LR on sampled data
#slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(np.log10(x), np.log10(random_zipf[a:b] ))
#print('LR slope for sampled data =', slope1)
#print('LR intercept for sampled data =', intercept1)

slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(np.log10(x), np.log10(random_zipf1[a:b] ))
print('LR_SF m for sampled data = ', slope2)
print('LR_SF q for sampled data = ', intercept2)
slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(np.log10(x), np.log10(random_zipf2[a:b] ))
print('LR_PMF m for sampled data = ', slope3)
print('LR_PMF q for sampled data = ', intercept3)


#Plotting sampled data and LR on sampled data, first sample
#y2 = x1**slope1 * 10**intercept1
#ax1.plot(x1,y2 , linewidth=1, color='purple', label = "LR on Zipf's sampled data")

#Plotting sampled data and LR on sampled data, second sample SF
#y3 = x1**slope2 * 10**intercept2
#ax1.plot(x1,y3 , linewidth=1, color='blue', label = "LR on Zipf's sampled data 2, SF")

#Plotting sampled data and LR on sampled data, second sample PMF
y4 = (np.arange(1, word+1))**slope3 * 10**intercept3
ax1.plot(np.arange(1, word+1),y4 , linewidth=1, color='green', label = "LR on Zipf's sampled data 2, PMF", zorder = 110)


#Setting log-log scale
ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.legend(loc='best')


ax1.set_title('Real data with their linear regressions plotted')

ax1.set_xlabel('Rank')
ax1.set_ylabel('Normalized frequency')

plt.show()


"""
Survival function: SF


#Survival function:
sf_y = np.zeros(word)
for i in range (0, word-1):
    sf_y[i] = norm_sum_y - np.sum(y[i+1:word])
    print(sf_y[i])
sf_y[word-1] =  norm_sum_y
sf_y = 1 - sf_y
#sf_y = sf_y/np.sum(sf_y)


#LR_SF on real data
x = np.arange(a, b)
sf_y = sf_y[a:b]
slope_sf, intercept_sf, r_value_sf, p_value_sf, std_err_sf = stats.linregress(np.log10(x), np.log10(sf_y))
print('LR_SF m for real data =', slope_sf)
print('LR_SF q for real data =', intercept_sf)


#plotting LR_SF for real data
c1_sf = 10**intercept_sf
x1_sf = np.arange(1, word+1)
y1_sf = (x1_sf**slope_sf) * c1_sf
#print(np.sum(y1))
#ax1.plot(x1_sf, y1_sf, label="LR_SF on real data", linewidth=1, color='xkcd:yellow brown')

"""