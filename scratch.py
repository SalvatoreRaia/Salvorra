import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import pickle
from scipy import optimize
from scipy import stats

def fa(x, a, n):
    return ((a-1)/(n-1)-(a-1)/(n-1)*(1/x))**(1/(1-a))*x**(1/(1-a))
def df1(x, a, n):
    return (((a - 1) * (x - 1) / n - 1) ** (a / (1 - a))) / (1 - n)
def df2(x,a,N):
    return (a* (((a - 1)*(x - 1))/(N - 1))**(1/(1 - a) - 2))/(N - 1)**2
# LINEAR REGRESSION (LR)
# Plotting function for linear regression

def lr(a, b, y, sample, label_lr = "Normalized Linear regression ", label_sample = "Zipf random sample "):
    # a = starting point, b = end point, y = da,
    # if sample == 1 sampled data will be plotted with same slope as the one found in LR.
    # x values fo LR are defined in [a+1,b]

    x_LR = np.arange(1+a , b+1)
    # Why 1+a : y[0] is defined in x = 1, so if we choose a = 0 then x = 1.
    # Why b+1: np.arange exclude b, so if I need b included, it's necessary ro put b+1 inside np.arange
    y_LR = y[a:b]

    #scipy LR
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x_LR), np.log10(y_LR))

    print('LR m for real data =', slope)
    print('LR q for real data =', intercept)

    # LR values, eq.3 report
    c = 10 ** intercept
    y_LR = (np.arange(a, a+len(y_LR)) ** round(slope, 3)) * c

    #plotting linear regression
    ax.plot(np.arange(a, a+len(y_LR)), y_LR/np.sum(y_LR), label= label_lr , linewidth=1.3, zorder=100, color ='orange')

    #plotting x_min Big Red Dot
    ax.scatter(a, y[a], s=30, color='red', label = '$x_{min}=$'+str(a)+ ", $\\alpha$ = " + str(round(slope, 3)) )

    if sample == 1:
        # Sampled data
        random_zipf1 = stats.zipf.pmf(np.arange(1,b+1), a= -slope)
        # Sample sorted
        random_zipf1 =  - np.sort(-random_zipf1)
        random_zipf2 = random_zipf1[a:b+1]/np.sum(random_zipf1[a:b+1])
        print(np.sum(random_zipf1), len(random_zipf1), len(random_zipf2))
        # Plotting
        ax.scatter(np.arange(1,len(random_zipf2)+1), random_zipf2, s=0.4, label=label_sample, color='green' )

    return y_LR[1:a+1]

x = np.arange(0,10000, 0.5)
y = 1/(60+x)**2

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Axis settings
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.grid(color='xkcd:light blue', linestyle='--', which='major', zorder=1, alpha=0.85)
ax.grid(color='xkcd:light blue', linestyle='--', which='minor', zorder=1, alpha=0.35)
ax.set_xticks(np.arange(0, 0.7, 0.05), minor=True)
ax.set_yticks(np.arange(0, 0.4, 0.05), minor=True)
# Setting log-log scale
ax.set_yscale('log')
ax.set_xscale('log')



#Which counties you want to be plotted? Also numpy array are ok.
#Chose a number between 1 and 3075
which_county = [3000]#, 69, 420, 666, 3000] #[1, 12, 75, 1356, 3000]
xmin = 22

#One of the row files to extract some usefull informations
Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header=None,
                    low_memory=False)


# total no. of counties
TOT_counties = len(Adata[1][:]) - 1

# List of total number of words per county, the order is the same as row data files
tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,
                       index_col=[0])

# Empty array to be filled with Total number of words of the i-th county
word = np.zeros(3074)

# CALLING FIGURE
# fig = plt.figure()
# ax1 = fig.add_subplot(1, 1, 1)

#IMPORTING DATA
with open(r'C:\Users\Salvo\Desktop\y.pkl', 'rb') as f:
    y = pickle.load(f)

print(y[0])

#CALCULATIONG CDF
cdf_y = [[]]
cdf_pd = [[]]
for i in range (3076):
    cdf_y.append([])
    cdf_pd.append([])

for c in which_county :

    county = Adata[1][c]

    dati = pd.read_csv("C:\\Users\Salvo\Desktop\IFISC\data\\row_frequencies\\frequencies_" + str(county) + "_(" + str(
        c) + ")" + ".txt", header=None).T
    #print(dat)
    dati = dati.to_numpy()[0]
    #print(dat, dat[0])

    dat = np.sort(np.sort(dati))
    unique, counts = np.unique(dat, return_counts=True)

    word[c] = int(tot_word.iloc[c, 1])
    # b is the total number of word for that county
    b = int(word[c])
    partial_sum = 0

    a = 1.49



    print('UNIQUE', len(unique))
    print('COUNTS', len(counts))
    print('y',len(y[c-1]))
    #print('COUNTS',counts)
    for j in range(b):
        partial_sum = partial_sum + (y[c-1])[j]
        cdf_y[c].append(partial_sum)
    b = len(dati)
    datipercdf = -np.sort(-dati)
    datipercdf = datipercdf/np.sum(datipercdf)
    print('SOMMA A 1 STA MERDA?', np.sum(datipercdf))
    partial_sum=0
    for j in range(b):
        partial_sum = partial_sum + datipercdf[j]
        cdf_pd[c].append(partial_sum)

        # scipy LR
    xmax=300

    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(unique[:xmax]), np.log10(counts[:xmax]))

    print('LR m for real data =', slope)
    print('LR q for real data =', intercept)

    # LR values, eq.3 report
    cc = 10 ** intercept
    y_LR = ((unique+100) ** round(slope, 3)) * cc

    counts = counts/np.sum(counts)
    #print('COUNTS NORMALIZED', counts, counts[0])
    xxx = np.arange(1,len(y[c-1])+1)
    ax.scatter(xxx, -np.sort(-dati/np.sum(dati)), s = 1.3, label =str(county)+' R-F')
    #ax.plot(-np.sort(-dati/np.max(dati)), -np.sort(-dati/np.max(dati)- (50+xxx )**-1.2), linewidth = 1.3, label =str(county)+' PDF')
    ax.scatter(unique, counts, s = 0.8, label = str(county+' PDF'))
    #ax.scatter(unique,(counts/np.max(counts)-unique**-1.5), s = 0.8, label = str(county+' R-F'))
    print('------------------------------\n a, b', a, b)

    xx = np.arange(len(xxx))
    yy = fa(xx, a, b)

    yyy = np.ones(len(xxx))
    for i in range(1,len(yy)):
        yy[i] = fa(xx[i],a,b)
    print(np.argmax(yy))
    #   yyy[i] = fa(xxx[i],a,b) + df1(xxx[i],a,b)*(xxx[i-1]-xxx[i]) + 0.5*df2(xxx[i],a,b)*(xxx[i-1]-xxx[i])**2
    #ax.scatter(xxx[::-1], yyy[::-1] / np.sum(yyy[::-1]), label='new eq', color='green', s = 0.8)
    ax.scatter(xx[1::], yy[1::]/np.sum(yy[1::]), label='new eq. true', color='purple', s = 0.8)


    # plotting linear regression
    ax.plot(unique, y_LR/np.sum(y_LR) , label='lr', linewidth=1.3, zorder=100, color='blue')
    print(len(unique[:xmax]))

ax.legend(loc='best')



figa = plt.figure()
ax1 = figa.add_subplot(1, 1, 1)


# Axis settings
ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.grid(color='xkcd:light blue', linestyle='--', which='major', zorder=1, alpha=0.85)
ax1.grid(color='xkcd:light blue', linestyle='--', which='minor', zorder=1, alpha=0.35)
ax1.set_xticks(np.arange(0, 0.7, 0.05), minor=True)
ax1.set_yticks(np.arange(0, 0.4, 0.05), minor=True)
# Setting log-log scale
#ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.legend(loc='best')

ax1.scatter(unique, counts - y_LR+1, label=str(county), color="orange", s =0.8)
ax1.plot([unique[0],unique[-1]], [1,1], color='blue')

#ax.plot(xxx, xxx**-1.2, label='suca')


 #LINEAR REGRESSION plot


plt.show()



