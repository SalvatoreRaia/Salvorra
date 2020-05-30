import numpy as np
import pandas as pd

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)
print(Adata)

#total no. of counties
N_counties = len(Adata[1][:])
print(N_counties)

#Adata = None

Alfabeto = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
dat = pd.DataFrame()

word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,
                   index_col=[0])

for i in range(len(Alfabeto)):
    rowfile = "C:\\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_" + Alfabeto[i] + "_PBW.txt"
    #print(rowfile)
    data = pd.read_csv(rowfile, sep=",", header=None, low_memory=False)
    data = data.drop(columns=[0, 2, 3])
    #print(data)

    for c in range (1, N_counties):
        county = Adata[1][c]
        #print(c)
        #print(county)

        dat = pd.concat([data.iloc[0, :], data.iloc[c, :]], axis=1)
        dat = dat.drop(index=[1])
        #print(dat)


        #FREQUENCIES
        freq_np = np.asarray(dat.iloc[:, 1].to_numpy(), dtype=float)
        a = np.flatnonzero(freq_np)
        freq_np = freq_np[a]
        #print(freq_np)

        f = open("C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\\row_frequencies\\frequencies_" + str(
            county) + "_(" + str(c) + ")" + ".txt", "a+")
        #f.write('\x08')
        #f.write(str(len(freq_np)) + '\n')

        for k in range(len(freq_np)):
            f.write(str(freq_np[k]) + '\n')
        f.close()



        #WORDS
        word_list = np.asarray(dat.iloc[:, 0].to_numpy())
        word_list = word_list[a]
        #print(word_list)

        w = open("C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\\list_of_words\\words_" +str(county)+"_("+ str(c)+")"+ ".txt", "a+")
        #w.write(str(len(word_list))+'\n')

        for k in range (len(word_list)):
            w.write(str(word_list[k])+'\n')
        w.close()

        print(len(Alfabeto)-i, N_counties-c)

    print(len(Alfabeto)-i)

print('fattoh')



