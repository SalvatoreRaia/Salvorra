import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)
prova = pd.DataFrame()


#total no. of counties
N_counties = len(Adata[1][:]) -1

#choose the county by a number between [1,N_counties]
c = 1234
county = list()
county.append( str(Adata[1][c]))
print(county)

N_counties = len(Adata[1][:]) -1


#chooose the number of words you need
word = 20


#how many counties you need?
hm_counties = 3

# GENERATE N INDEX THAT IDENTIFIES N COUNTIES.
# N_counties+2 because the idx of last county is N_counties+1, and rand int is [min, max)
for i in range(1, hm_counties):
    a = np.random.randint(1, N_counties + 2)
    c = np.append(c, a)
    print(c)
    county.append(str(Adata[1][c[i]]))
    print(county)




Adata = None

Alfabeto = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
dat = pd.DataFrame()
data2 = pd.DataFrame()

for k in range (0, hm_counties):
    for i in range(len(Alfabeto)):
        str = "C:\\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_" + Alfabeto[i] + "_PBW.txt"
        print(str)
        data = pd.read_csv(str, sep=",", header=None, low_memory=False)
        data = data.drop(columns=[0, 2, 3])
        print(data)

        # TAKINNG ALL FREQUENCIES
        data2 = data.iloc[c[k], :]
        # TAKING ALL WORDS
        data = data.iloc[0, :]

        # TRANSPOSING TO COLUMNS
        data = data.T
        data2 = data2.T

        # MERGING THE 2 DATAFRAMES
        data = pd.concat([data, data2], axis=1)

        data = data.drop(index=[1])
        data[c[k]] = pd.to_numeric(data[c[k]])
        print(data)

        # SORTING, RESET INDEX, TAKING FIRST N WORDS
        data = data.sort_values(by=[c[k]], ascending=False)
        data = data.reset_index()
        data = data.iloc[:word, :]

        # APPENDING TO THE FINAL DATAFRAME
        dat = dat.append(data, ignore_index=True)
        print(dat)

    # TAKING FIRST N WORDS OF THE FINAL DATAFRAME
    dat = dat.sort_values(by=[c[k]], ascending=False)
    dat = dat.reset_index()
    dat = dat.iloc[:word, :]
    print(dat)

    # GRAPHIC
    x = np.arange(1, word + 1)
    y = dat.iloc[:, 3]
    parole = dat.iloc[:, 2]
    print(x, y)

    fig = plt.figure()
    fig.set_dpi(100)
    ax1 = fig.add_subplot(1, 1, 1)

    xy = (x + 0.3, y + 0.3)
    ax1.scatter(x, y, label=county[k], color='red')
    # ax1.annotate(parole, xy, fontsize=9)
    ax1.legend(loc='upper right')
    ax1.grid(color='xkcd:cherry', linestyle='--', which='major', zorder=1, alpha=0.75)
    ax1.set_xticks(np.arange(0, word + 1), minor=False)
    # ax1.set_yticks(np.arange(0, 1, 0.1), minor=False)
    for i, parole in enumerate(parole):
        ax1.annotate(parole, (x[i] - 0.2, y[i] + 1.2))

    data2 = pd.DataFrame()
    data = pd.DataFrame()
    dat = pd.DataFrame()

    print(dat, data)

plt.show()



