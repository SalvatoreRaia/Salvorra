import numpy as np
import pandas as pd
import geopy.distance as gpdis

Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header=None,
                    low_memory=False)
#total no. of counties
TOT_counties = len(Adata[1][:]) -1

print(Adata)

data = Adata.drop(columns=np.arange(4,3801))
data = data.drop(columns=[0,1])
print(data)

long = np.asarray(data.iloc[1:, 0].to_numpy(), dtype=float)
lati = np.asarray(data.iloc[1:, 1].to_numpy(), dtype=float)

print('long: ',long)
print('lat: ',lati)

n_counties = TOT_counties
S = 3000
grid = np.arange(1, n_counties)
D = np.zeros((len(grid), len(grid)))
combinations = len(grid)*(len(grid)-1)*0.5

"""
Probably the filling doesn't include the last county.
I added it manually afterward. This should be fixed but the distance calulation for
a single couple of counties is ok
"""
#Filling triangular matrix
for i in range (0, len(grid)):
    for j in range (i+1, len(grid)):
        D[i][j] = gpdis.geodesic( (lati[i],long[i]) , (lati[j], long[j]) ).km
        #np.savez(r"C:\Users\Salvo\Desktop\IFISC\data\all_counties\ks\D_pvalue_(" + str(i) + '-' + str(j) + ")", counties=str(i)+','+str(j), D=D[i][j])

    print(len(grid) - i )

np.savez(r'C:\Users\Salvo\Desktop\IFISC\data\all_counties\g_distances\run\g_D_'+str(int(combinations))+'_combinations', D=D,  counties = grid)

print('fattoh')