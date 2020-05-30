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



#One of the row files to extract some usefull informations
Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header=None,
                    low_memory=False)
print(Adata)

#cutting out some columns
Adata = Adata.drop(columns=np.arange(4,3801))
Adata = Adata.drop(columns=[0,1])
print(Adata)

#taking latitude and longitude of all the counties, it is like an avarage point within a county
long = np.asarray(Adata.iloc[1:-1, 0].to_numpy(), dtype=float)
lati = np.asarray(Adata.iloc[1:-1, 1].to_numpy(), dtype=float)

#KS
#Importing KS
combinations = 4723201
myks_file = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\all_counties\ks\run\myD_value_'+str(int(combinations))+'_combinations.npz')

#Extracting arrays from file
myD = myks_file['D']
ks_grid = myks_file['counties']

#This is one of the shuffled ks distance matrix
#myks_file = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\New_Dist\N_ks_D_7.npz')
#myD = myks_file['N_ks_D']
#ks_grid = myks_file['N_grid']

print(myD, len(myD[0][:]))


#JS
combinations = 4723201
#Importing JS
js_distance = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\all_counties\js\run\jsD_value_'+str(int(combinations))+'_combinations.npz')
print(js_distance.files)

#Extracting arrays from file
js_D = js_distance['jd_D']
js_grid = js_distance['counties']

#js_distance = np.load(r'C:\Users\Salvo\Desktop\IFISC\data\New_Dist\N_js_D_2.npz')
#js_D = js_distance['N_js_D']
#js_grid = js_distance['N_grid']

# List of total number of words per county
tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,index_col=[0])
print(tot_word)
word = np.asarray(tot_word.iloc[1:-1, 1].to_numpy(), dtype=float)
print(word)

#Plotting a scatter plot for each number of clusters
which_clusters = [3,5,7,9]


#"""
#KS
for i in which_clusters:

    #calling the object
    ks_cluster = AgglomerativeClustering(n_clusters=i, affinity='euclidean', linkage='ward')
    a = ks_cluster.fit_predict(myD)
    print(a)
    #a has in position i the county with index i.
    #as value in position i , it has the number of the cluster that the i-th county belongs to (example, 3 clusters, a =[1,1,1,1,0,0,0,2,1,1]

    #Here it is extracting the how many counties are in each cluster, (example, 3 clusters, bin_count = [0, 3, 6, 1])
    unique, bin_count_temp = np.unique(a, return_counts=True)
    #with concatenate I put a 0 as a first element of bin count
    bin_count = np.concatenate(([0], bin_count_temp))

    print(bin_count, len(bin_count))

    #b are the counties indexes
    b = np.arange(len(a))

    #Here I create a pandas dataframe to easily sort counties, latitude and longitude arrays, following the clusters sorting.
    ok = pd.DataFrame({'cluster': a,'counties':b, 'latitudine':lati, 'longitudine':long})
    ok = ok.sort_values(by=['clusters'])
    print(ok)

    #I will need his after
    sorted_count = np.sort(bin_count)

    #to dataframe to array
    long = np.asarray(ok.iloc[:, 3].to_numpy(), dtype=float)
    lati = np.asarray(ok.iloc[:, 2].to_numpy(), dtype=float)


    #Here I create a figure for each number of cluster
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    # now i want to draw a scatter plot x = longitude, y = latitude with different clusters in different colors
    for k in range (len(bin_count)-1): #


        #from here
        da = np.sum(bin_count[:k + 1])
        #to here
        aa = da + bin_count[k + 1]

    #plot those points with the same colors as they belongs to the same cluster
    #But just the 4 biggest clusters, all the others in grey

        #if clusters are less than 4
        if len(sorted_count) < 4:
            ax.scatter(long[da:aa], lati[da:aa], s=3.5, label=str(aa - da))


        else:

            #if the number of counties in this k-th cluster are less than the fourth biggest cluster
            fourth = sorted_count[-4]
            if (aa - da) < (fourth):
                #color it in grey
                ax.scatter(long[da:aa], lati[da:aa], s=0.6, color='grey', label=str(aa - da))
            else:
                ax.scatter(long[da:aa], lati[da:aa], s=3.5, label=str(aa - da))
        ax.set_title('KS ' + str(i) + ' clusters')

    print('KS con:', i)
    #writing how  many grey counties in the legend of the graph
    if len(sorted_count) > 4:
        ax.scatter(-128, 25, s=0.001, label='Grey counties: ' + str(np.sum(sorted_count[:-4])))
        print('Grey counties', np.sum(sorted_count[:-4]))
    ax.legend(loc='best')
#"""

#"""
#JS
for i in which_clusters:
    ks_cluster = AgglomerativeClustering(n_clusters=i, affinity='euclidean', linkage='ward')
    a = ks_cluster.fit_predict(js_D)
    unique, bin_count_temp = np.unique(a, return_counts=True)
    bin_count = np.concatenate(([0], bin_count_temp))

    print(bin_count, len(bin_count))
    b = np.arange(len(a))
    ok = pd.DataFrame({'clusters': a,'counties':b, 'latitudine':lati, 'longitudine':long, 'parole': word})
    ok = ok.sort_values(by=['clusters'])
    print(ok)

    sorted_count = np.sort(bin_count)

    long = np.asarray(ok.iloc[:, 3].to_numpy(), dtype=float)
    lati = np.asarray(ok.iloc[:, 2].to_numpy(), dtype=float)
    counties = np.asarray(ok.iloc[:, 1].to_numpy(), dtype=int)
    word = np.asarray(ok.iloc[:, 4].to_numpy(), dtype=int)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)


    for k in range (len(bin_count)-1):
        da = np.sum(bin_count[:k+1])
        aa = da + bin_count[k+1]
        print(da, aa)

        #ax1.plot(counties[da:aa],word[da:aa] , linewidth = 0.5, label = 'da '+str(da)+' a '+str(aa))

        if len(sorted_count) < 4:
            ax.scatter(long[da:aa], lati[da:aa], s=3.5, label=str(aa - da))
        else:
            # ax.scatter(-128, 25, s = 0.001, label = 'Contee in grigio: '+ str(np.sum(sorted_count[:-3])))
            third = sorted_count[-4]
            if (aa - da) < (third):
                ax.scatter(long[da:aa], lati[da:aa], s=0.6, color='grey', label=str(aa - da))
            else:
                ax.scatter(long[da:aa], lati[da:aa], s=3.5, label=str(aa - da))
        ax.set_title('JS ' + str(i) + ' clusters')

    print('JS con:', i)
    if len(sorted_count) > 4:
        ax.scatter(-128, 25, s=0.001, label='Contee in grigio: ' + str(np.sum(sorted_count[:-4])))
        print('Contee in grigio', np.sum(sorted_count[:-4]))
    ax.legend(loc='best')


plt.show()
