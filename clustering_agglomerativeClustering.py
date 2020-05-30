from sklearn.cluster import AgglomerativeClustering
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#########################################################
"""
https://scikit-learn.org/stable/modules/clustering.html
documentation for clustering algorithm
"""
#########################################################


# One of the row files to extract some usefull informations
Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header=None,
                    low_memory=False)
print(Adata)

# cutting out some columns
Adata = Adata.drop(columns=np.arange(4, 3801))
Adata = Adata.drop(columns=[0, 1])
print(Adata)

# taking latitude and longitude of all the counties, it is like an avarage point within a county
long = np.asarray(Adata.iloc[1:, 0].to_numpy(), dtype=float)
lati = np.asarray(Adata.iloc[1:, 1].to_numpy(), dtype=float)

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

# This is one of the shuffled ks distance matrix
# myks_file = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\New_Dist\N_ks_D_7.npz')
# myD = myks_file['N_ks_D']
# ks_grid = myks_file['N_grid']

# shuffled matrix
# js_distance = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\New_Dist\N_js_D_2.npz')
# js_D = js_distance['N_js_D']
# js_grid = js_distance['N_grid']

# List of total number of words per county
tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,
                       index_col=[0])
print(tot_word)
word = np.asarray(tot_word.iloc[1:, 1].to_numpy(), dtype=float)
print(word)

# Plotting a scatter plot for each number of cluster
which_clusters = [10, 12, 20, 40]  # , 7, 9]

# First n biggest clusters will be colorized
n = 50

#################################
# KS
#################################

for i in which_clusters:

    # calling the object
    ks_cluster = AgglomerativeClustering(n_clusters=i, affinity='euclidean', linkage='ward')
    a = ks_cluster.fit_predict(myD)
    print(len(a), a, a[12])
    # a has in position i the county with index i. as value in position i , it has the number of the cluster that the
    # i-th county belongs to (example, 3 clusters, a =[1,1,1,1,0,0,0,2,1,1,...,2], a[3] = 1, this means that the fourth
    # county is ine the cluster 1

    # Here it is extracting the how many counties are in each cluster.
    # (example, 3 clusters, bin_count=[0, 1478, 1592, 5]) 1478 counties in the 0th county, 1592 in the 2th, 5 in the 3th
    unique, bin_count_temp = np.unique(a, return_counts=True)
    # with concatenate I put a 0 as a first element of bin count
    bin_count = np.concatenate(([0], bin_count_temp))
    print(bin_count, len(bin_count))

    # Sorting it
    sorted_count = np.sort(bin_count)  # I will need this later

    # b are the counties indexes
    b = np.arange(len(a))

    # Here I create a pandas dataframe to easily sort counties, latitude and longitude arrays, following the clusters
    # sorting.
    ok = pd.DataFrame({'clusters': a, 'counties': b, 'latitudine': lati, 'longitudine': long})
    ok = ok.sort_values(by=['clusters'])
    print(ok)

    # Dataframe to array
    long = np.asarray(ok.iloc[:, 3].to_numpy(), dtype=float)
    lati = np.asarray(ok.iloc[:, 2].to_numpy(), dtype=float)

    # Here I create a figure for each number of cluster
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # now i want to draw a scatter plot x = longitude, y = latitude with different clusters in different colors
    for k in range(len(bin_count) - 1):  #

        # from here
        da = np.sum(bin_count[:k + 1])
        # to here
        aa = da + bin_count[k + 1]

        # plot those points with the same colors as they belongs to the same cluster
        # But just the 4 biggest clusters, all the other smaller clusters in grey

        # if clusters are less than 4
        if len(sorted_count) < n:
            ax.scatter(long[da:aa], lati[da:aa], s=15, label=str(aa - da))

        else:

            # if the number of counties in this k-th cluster are less than the n-th biggest cluster
            nth = sorted_count[-n]
            if (aa - da) < nth:
                # color it in grey
                ax.scatter(long[da:aa], lati[da:aa], s=9, color='grey', label=str(aa - da))
            else:
                ax.scatter(long[da:aa], lati[da:aa], s=15, label=str(aa - da))
        ax.set_title('KS ' + str(i) + ' clusters')

    print('KS con:', i)
    # writing how  many grey counties in the legend of the graph
    if len(sorted_count) > n:
        ax.scatter(-128, 25, s=0.001, label='Grey counties: ' + str(np.sum(sorted_count[:-n])))
        print('Grey counties', np.sum(sorted_count[:-n]))
    ax.legend(loc='best')


#################################
# JS
#################################

for i in [2,10]:#which_clusters:
    ks_cluster = AgglomerativeClustering(n_clusters=i, affinity='euclidean', linkage='ward')
    a = ks_cluster.fit_predict(js_D)
    unique, bin_count_temp = np.unique(a, return_counts=True)
    bin_count = np.concatenate(([0], bin_count_temp))

    print(bin_count, len(bin_count))
    b = np.arange(len(a))
    ok = pd.DataFrame({'clusters': a, 'counties': b, 'latitudine': lati, 'longitudine': long, 'parole': word})
    ok = ok.sort_values(by=['clusters'])
    print(ok)

    sorted_count = np.sort(bin_count)

    long = np.asarray(ok.iloc[:, 3].to_numpy(), dtype=float)
    lati = np.asarray(ok.iloc[:, 2].to_numpy(), dtype=float)
    counties = np.asarray(ok.iloc[:, 1].to_numpy(), dtype=int)
    word = np.asarray(ok.iloc[:, 4].to_numpy(), dtype=int)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    for k in range(len(bin_count) - 1):
        da = np.sum(bin_count[:k + 1])
        aa = da + bin_count[k + 1]
        print(da, aa)

        # ax1.plot(counties[da:aa],word[da:aa] , linewidth = 0.5, label = 'da '+str(da)+' a '+str(aa))

        if len(sorted_count) < n:
            ax.scatter(long[da:aa], lati[da:aa], s=15, label=str(aa - da), marker='s')
        else:
            # ax.scatter(-128, 25, s = 0.001, label = 'Contee in grigio: '+ str(np.sum(sorted_count[:-3])))
            nth = sorted_count[-n]
            if (aa - da) < nth:
                ax.scatter(long[da:aa], lati[da:aa], s=9, color='grey', label=str(aa - da), marker='s')
            else:
                ax.scatter(long[da:aa], lati[da:aa], s=15, label=str(aa - da), marker='s')
        ax.set_title('JS ' + str(i) + ' clusters')

    print('JS con:', i)

    # If there are more than 5 clusters, a label with the total number of grey clusters will be plotted
    if len(sorted_count) > n:
        ax.scatter(-128, 25, s=0.001, label='Contee in grigio: ' + str(np.sum(sorted_count[:-n])))
        print('Contee in grigio', np.sum(sorted_count[:-n]))
    ax.legend(loc='best')

plt.show()
