import numpy as np
import pandas as pd


def gen_county(c, word = None):
    Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)

    #total no. of counties
    N_counties = len(Adata[1][:]) -1
    print(N_counties)

    #choose the county by a number between [1,N_counties]
    #c = np.random.choice(np.arange(1, N_counties + 1 ))
    #c = 273
    county = Adata[1][c]
    print(c)
    print(county)

    #chooose the number of words you need; default is the maximum
    if word is None :
        word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False, index_col=[0])
        word = int(word.iloc[c, 1]) #+1 is due to the indexing ':word' that don't take last word
        #word = 3 #custom number
        print(word)
    Adata = None

    Alfabeto = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
    dat = pd.DataFrame()

    for i in range(len(Alfabeto)):
        rowfile = "C:\\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_" + Alfabeto[i] + "_PBW.txt"
        print(rowfile)
        data = pd.read_csv(rowfile, sep=",", header=None, low_memory=False)
        data = data.drop(columns=[0, 2, 3])
        print(data)



        # TAKINNG ALL FREQUENCIES
        data2 = data.iloc[c, :]
        # TAKING ALL WORDS
        data = data.iloc[0, :]

        # TRANSPOSING TO COLUMNS
        data = data.T
        data2 = data2.T

        # MERGING THE 2 DATAFRAMES
        data = pd.concat([data, data2], axis=1)

        data = data.drop(index=[1])
        data[c] = pd.to_numeric(data[c])
        print(data)

        # SORTING, RESET INDEX, TAKING FIRST N WORDS
        data = data.sort_values(by=[c], ascending=False)
        data = data.reset_index()
        data = data.iloc[:word, :]

        # APPENDING TO THE FINAL DATAFRAME
        dat = dat.append(data, ignore_index=True)
        print(dat)

    # TAKING FIRST N WORDS OF THE FINAL DATAFRAME
    dat = dat.sort_values(by=[c], ascending=False)
    dat = dat.reset_index()
    dat = dat.iloc[:word+26, :]
    print(dat)

    # SAVING DATA ON FILE
    z = np.zeros(word+26)
    z[0] = int(c)

    np.savez("C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\\row_frequencies\\frequencies_" + str(
        county) + "_(" + str(c) + ")", freq = dat.iloc[:, 3])

    rfcw = pd.DataFrame({'words': dat.iloc[:, 2],'frequency':dat.iloc[:, 3], 'county_(idx_no.)': z})
    print(rfcw)
    filename = r"C:\Users\Salvo\Desktop\IFISC\data\r_f_county_w\freq_ranking_"+county+'(' + str(c) +')_'+str(word)+"w.csv"
    rfcw.to_csv(filename)

    print('fattoh')

#choose the county by a number between [1,N_counties(3075)]
c = np.random.choice(np.arange(1, 3075 + 1 ))
c = 2604

gen_county(c)
print(c)




