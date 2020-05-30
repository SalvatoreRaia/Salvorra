import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

zipf_sample= np.random.zipf(1.1, size=(11000))

zipf_sample = - np.sort(-zipf_sample)
print(zipf_sample)

x = np.arange(1,11001)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.scatter(x,zipf_sample, s=1 )

ax.set_yscale('log')
ax.set_xscale('log')

slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x), np.log10(zipf_sample[a:(b + a - i)]))











'''

a = 200
b = 8000

sloppete = np.zeros(b-a+1)
interceptete = np.zeros(b-a+1)
std_errete = np.zeros(b-a+1)

for i in range(a, b, 1):
    x = np.arange(a, b+a-i)

    print(x, zipf_sample)
    print(b - i)
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x), np.log10(zipf_sample[a:(b+a-i)]))

    sloppete[b-i] = slope
    interceptete[b-i] = intercept
    std_errete[b-i] = std_err

    #x1 = np.arange(1, word + 1, 0.1)
    #y1 = 10**intercept*x1**(slope)
    #ax1.plot(x1,y1)

    #print(slope,10**(intercept), std_err)
    #print(parole)

x = np.arange(a, b+1)
#x = np.arange(a, b, 1)
ax2 = fig.add_subplot(1, 1, 1)
ax2.scatter(x, sloppete, s = 0.5 )
#ax2.scatter(x, std_errete,s = 0.5)
#filename = r"C:\Users\Salvo\Desktop\IFISC\data\graphics\freq_ranking_"+county+'(' + str(c) +')_'+str(word)+"w.eps"
#fig.savefig('C:\\Users\Salvo\Desktop\IFISC\data\graphics', format='eps')

plt.show()


'''






'''
word = 30
c = 1249
county = 'bauabu'

# GRAPHIC
x = np.arange(1, word + 1)
y = x**-1.1
print(x, y)

fig = plt.figure()
fig.set_dpi(100)
ax1 = fig.add_subplot(1, 1, 1)

ax1.scatter(x, y, label=county, color='red')
ax1.legend(loc='upper right')
ax1.grid(color='xkcd:cherry', linestyle='--', which='major', zorder=1, alpha=0.75)
ax1.set_xticks(np.arange(0, word + 1), minor=False)
#for i, parole in enumerate(parole):
 #   ax1.annotate(parole, (x[i] - 0.2, y[i] + 1.2))

ax1.set_yscale('log')
ax1.set_xscale('log')

#Linear regression
slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x),np.log10(y))


# creating a raking, frequency, county array
z = np.zeros(len(x))
print(c)
z[0] = c


parole = pd.DataFrame({'words':["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T"]})
rfc = pd.DataFrame({'x_values':x[:], 'y_values':y[:], 'county_number':z[:]})
rfcw = pd.concat([rfc, parole], axis=1)
print(rfcw)
filename = r"C:\\Users\Salvo\Desktop\IFISC\data\\r_f_county_w\\freq_ranking_"+county+str(word)+"w.csv"
rfcw.to_csv(filename)


provaimporta = pd.read_csv(filename, sep=",", header= None, low_memory=False)
provaimporta = provaimporta.drop(columns=[0])
provaimporta = provaimporta.reset_index()
print(provaimporta)


#print(slope,10**(intercept))
#print(rfc)


plt.show()

'''
