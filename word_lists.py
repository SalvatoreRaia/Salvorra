import numpy as np
import pandas as pd

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)

#total no. of counties
N_counties = len(Adata[1][:]) -1
print(N_counties)

#Adata = None

Alfabeto = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
dat = pd.DataFrame()

for i in range(len(Alfabeto)):
    rowfile = "C:\\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_" + Alfabeto[i] + "_PBW.txt"
    print(rowfile)
    data = pd.read_csv(rowfile, sep=",", header=None, low_memory=False)
    data = data.drop(columns=[0, 2, 3])
    print(data)

    data = data.iloc[0][1:]
    data.reset_index()
    print(data)

    # SAVING DATA ON FILE
    #z = np.zeros(len(word_list))
    #rfcw = pd.DataFrame({'words': dat.iloc[:, 2], 'frequency': dat.iloc[:, 3], 'county_(idx_no.)': z})
    #print(rfcw)
    filename = "C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\\" + Alfabeto[i] + "_word_list.csv"
    data.to_csv(filename)

    print('fattoh')
'''
    for c in range (10):
        # choose the county by a number between [1,N_counties]
        # c = np.random.choice(np.arange(1, N_counties + 1 ))
        # c = 273

        county = Adata[1][c]
        print(c)
        print(county)


        # chooose the number of words you need; default is the maximum
        word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,
                           index_col=[0])
        word = int(word.iloc[c, 1])  # +1 is due to the indexing ':word' that don't take last word
        # #word = 150 custom number
        print(word)

'''
