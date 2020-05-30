import numpy as np
import pandas as pd

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)
prova = pd.DataFrame()

#total no. of counties
N_counties = len(Adata[1][:]) -1
print(N_counties)

#counting (words) per county
county_list = Adata.iloc[:, 1]
counting_list = np.zeros(N_counties+1)

Alfabeto = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]

for i in range(len(Alfabeto)):
    rowfile = "C:\\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_" + Alfabeto[i] + "_PBW.txt"
    print(rowfile)
    data = pd.read_csv(rowfile, sep=",", header=None, low_memory=False)
    data = data.drop(columns=[0, 2, 3])
    print(data)
    tot_parole = len(data.iloc[0, :]) - 1
    flagg = np.zeros(len(data.iloc[0, :]) - 1)
    #data = data.T
    #print(data)
    for k in range(1, N_counties + 1):
        print(N_counties + 1 - k)
        a = data.iloc[k, 1:].to_numpy()
        print(a)
        a = a.astype(np.float)
        a, = np.nonzero(a)
        a = len(a)
        print(a)
        counting_list[k] = counting_list[k] + a
        print(counting_list)
    print(counting_list)

cp_county = pd.DataFrame({'counting':counting_list})
wp_county = pd.concat([county_list, cp_county], axis=1)
print(wp_county)

#filename = r"C:\Users\Salvo\Desktop\IFISC\data\no_words_per_county.csv"
#wp_county.to_csv(filename)
