from multiprocessing import Pool
import pickle
import numpy as np
from itertools import product




def js_distance(p, q):
    # checking leght
    a = len(p) - len(q)
    if a < 0:
        p = np.append(p, np.zeros(-a))
    elif a > 0:
        q = np.append(q, np.zeros(a))

    # m = mean p,q
    m = (p + q) * 0.5

    kl_D1 = np.sum(p * np.log(p / m))
    kl_D2 = np.sum(q * np.log(q / m))

    # print('1:',kl_d)
    if np.isnan(kl_D1):
        kl_D1 = 0
    if np.isnan(kl_D2):
        kl_D2 = 0

    # js distance
    js_d = np.sqrt(0.5 * (kl_D1 + kl_D2))
    print(js_d)
    return js_d



#Filling triangular matrix
#for i in range (0, len(y[:])):
#    for j in range (i+1, len(y[:])):
#        iterabbile.append(())

if __name__ == '__main__':
    with open(r'C:\Users\Salvo\Desktop\y.pkl', 'rb') as f:
        y = pickle.load(f)

    print(y)
    print(len(y[0][:]))

    contee = np.arange(3075)
    with Pool(8) as p:
        results = p.map(js_distance, product(y, repeat=2))
        print(results)
