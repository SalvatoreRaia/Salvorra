import pandas as pd
import numpy as np

pop = pd.read_csv(r'C:\Users\Salvo\Desktop\co-est2018-alldata.csv', sep=",", header= None, low_memory=False, encoding = "latin 1")
print(pop)

#print(pop.iloc[:, 5])
#print(pop.iloc[:, 6])
#print(pop.iloc[:, 7])

pop = pd.DataFrame({str(pop.iloc[0, 5]
): pop.iloc[1:, 5], str(pop.iloc[0, 6]
): pop.iloc[1:, 6], str(pop.iloc[0, 7]
): pop.iloc[1:, 7]
})

print(pop)

filename = r"C:\Users\Salvo\Desktop\IFISC\data\population.csv"
pop.to_csv(filename)

print('fattoh')

popa = pd.read_csv(filename, sep=",", header= None, low_memory=False)
print(popa)

#COUNTIES FROM POP DATA
counties = list(sorted(set(popa.iloc[:, 2].tolist())))
print('counties',counties)
for i in range (len(counties)):
    counties[i] = counties[i].lower()
print('should be lower counties',counties)

#COUNTIRES
uni_countries = list(sorted(set(popa.iloc[:, 1].tolist())))
print('Uni countries',uni_countries)

#indices = [i for i, s in enumerate(counties) if 'County' in s]
#no_index = [i for i, s in enumerate(counties) if  not 'County' in s]
#print(no_index)
#print(len(no_index))
#print(indices)
#print(len(indices))
#print(len(no_index)+len(indices), len(counties))
#print(popa.iloc[:, 1].tolist())




#MY DATA

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header= None, low_memory=False)
print(Adata)

county_list = Adata.iloc[:,1].to_list()
lista1 = []
lista2 = []
print(county_list)
print('county list my data:',len(county_list))

indexes = []
indexess = []
lenghtcounties= []
lenghtcountiess= []

for i in range (1,len(county_list)):
    a = county_list[i]
    indexes.extend([pos for pos, char in enumerate(a) if char == ';' ])
    #indexess.extend([pos for pos, char in enumerate(a) if char == ':' ])
    lista2.append(county_list[i][  indexes[i-1]+1:])
    #lista1.append(county_list[i][ : indexes[i-1]])
    lenghtcountiess.append(len(lista2[i-1]))


print('lunghezza delle contee:',lenghtcountiess)
print(lista2)

#lista1 = sorted(set(lista1))
#lista2 = list(lista2)
#print(lista2)


'''
aa = ['a', 'b', 'c']
bb = ['c', 'b', 'd', 'a']

[a for a, b in zip(aa, bb) if a==b]
'''

popu_county = np.zeros(len(lista2))

for i in range (len(lista2)) :
    print(len(lista2)-i)
    for j in range (1, len(counties)):
        a = str(counties[j])
        #print(a)
        if str(lista2[i]) == str(a[:len(lista2[i])]) :
            popu_county[i] =  float(popa.iloc[j, 3])
        else :
            popu_county[i] = None

print(popu_county)
