import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from matplotlib import pyplot as plt

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header=None,
                    low_memory=False)
print(Adata)

Adata = Adata.drop(columns=np.arange(4, 3801))
Adata = Adata.drop(columns=[0, 1])
print(Adata)

long = np.asarray(Adata.iloc[1:, 0].to_numpy(), dtype=float)
lati = np.asarray(Adata.iloc[1:, 1].to_numpy(), dtype=float)
print(len(Adata))

print(long, len(long))
print(lati, len(lati))

# KS
# loading distance matrix
npzfile2 = np.load(r'C:\Users\Salvo\Desktop\Project\ks_distanceMatrix.npz')
print(npzfile2.files)
# assigning arrays in the file to new arrays
myD = npzfile2['ks_D']
# grid is the list of the counties in the order they appear in the ditsance matrix
ks_grid = npzfile2['counties']

# JS
js_distance = np.load(r'C:\Users\Salvo\Desktop\Project\js_distanceMatrix.npz')
print(js_distance.files)
js_D = js_distance['js_D']
js_grid = js_distance['counties']

# List of total number of words per county
tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,
                       index_col=[0])
print(tot_word)
word = np.asarray(tot_word.iloc[1:, 1].to_numpy(), dtype=float)
print(word)

# Plotting a scatter plot for each number of cluster
which_clusters = [10]  # [3, 5, 7, 9, 13, 17, 25, 41]

# First n biggest clusters will be colorized
n = 5

# '''
#################################
# KS
#################################

for i in which_clusters:
    ks_cluster = KMeans(n_clusters=i)
    a = ks_cluster.fit_predict(myD)
    unique, bin_count_temp = np.unique(a, return_counts=True)
    bin_count = np.concatenate(([0], bin_count_temp))

    sorted_count = np.sort(bin_count)
    # print('sorted count:',sorted_count)

    # print(bin_count, len(bin_count))
    b = np.arange(len(a))
    ok = pd.DataFrame({'clusters': a, 'counties': b, 'latitudine': lati, 'longitudine': long})
    ok = ok.sort_values(by=['clusters'])
    # print(ok)

    long = np.asarray(ok.iloc[:, 3].to_numpy(), dtype=float)
    lati = np.asarray(ok.iloc[:, 2].to_numpy(), dtype=float)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for k in range(len(bin_count) - 1):
        da = np.sum(bin_count[:k + 1])
        aa = da + bin_count[k + 1]
        print(da, aa)

        if len(sorted_count) < n:
            ax.scatter(long[da:aa], lati[da:aa], s=15, label=str(aa - da))
        else:
            nth = sorted_count[-n]
            if (aa - da) < nth:
                ax.scatter(long[da:aa], lati[da:aa], s=9, color='grey', zorder=100, label=str(aa - da))
            else:
                ax.scatter(long[da:aa], lati[da:aa], s=15, label=str(aa - da))
        ax.set_title('KS ' + str(i) + ' clusters')

    print('KS con: ', i)

    if len(sorted_count) > n:
        ax.scatter(-128, 25, s=0.001, label='Contee in grigio: ' + str(np.sum(sorted_count[:-n])))
        print('Contee in grigio', np.sum(sorted_count[:-n]))

    # Graphic Settings

    # Legend
    ax.legend(loc='best')

    # Axis grid
    ax.grid(color='xkcd:light blue', linestyle='--', which='major', zorder=1, alpha=0.85)
    ax.grid(color='xkcd:light blue', linestyle='--', which='minor', zorder=1, alpha=0.35)
    ax.set_xticks(np.arange(1, 70000, 100), minor=True)
    ax.set_yticks(np.arange(1, 70000, 100), minor=True)

    # Axis titles
    ax.set_xlabel('Ranking')
    ax.set_ylabel('Normalized frequency')

# '''
#################################
# JS
#################################

for i in which_clusters:
    ks_cluster = KMeans(n_clusters=i)
    a = ks_cluster.fit_predict(js_D)
    unique, bin_count_temp = np.unique(a, return_counts=True)
    bin_count = np.concatenate(([0], bin_count_temp))

    sorted_count = np.sort(bin_count)

    # print(bin_count, len(bin_count))
    b = np.arange(len(a))
    ok = pd.DataFrame({'clusters': a, 'counties': b, 'latitudine': lati, 'longitudine': long, 'parole': word})
    ok = ok.sort_values(by=['clusters'])
    # print(ok)

    long = np.asarray(ok.iloc[:, 3].to_numpy(), dtype=float)
    lati = np.asarray(ok.iloc[:, 2].to_numpy(), dtype=float)
    counties = np.asarray(ok.iloc[:, 1].to_numpy(), dtype=int)
    word = np.asarray(ok.iloc[:, 4].to_numpy(), dtype=int)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # fig1 = plt.figure()
    # ax1 = fig1.add_subplot(1, 1, 1)

    for k in range(len(bin_count) - 1):
        da = np.sum(bin_count[:k + 1])
        aa = da + bin_count[k + 1]
        # print(da, aa)

        # ax1.scatter(counties[da:aa], word[da:aa], s=0.5, label='da ' + str(da) + ' a ' + str(aa))

        if len(sorted_count) < n:
            ax.scatter(long[da:aa], lati[da:aa], s=15, label=str(aa - da))
        else:
            # ax.scatter(-128, 25, s=0.001, label='Contee in grigio: ' + str(np.sum(sorted_count[:-3])))
            nth = sorted_count[-n]
            if (aa - da) < nth:
                ax.scatter(long[da:aa], lati[da:aa], s=9, color='grey', zorder=100, label=str(aa - da))
            else:
                ax.scatter(long[da:aa], lati[da:aa], s=15, label=str(aa - da))

        ax.set_title('JS ' + str(i) + ' clusters')

    print('JS con:', i)

    if len(sorted_count) > n:
        ax.scatter(-128, 25, s=0.001, label='Contee in grigio: ' + str(np.sum(sorted_count[:-n])))
        print('Contee in grigio', np.sum(sorted_count[:-n]))

    ax.legend(loc='best')

plt.show()
