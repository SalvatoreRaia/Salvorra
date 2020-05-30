


# define a example function
if __name__ == '__main__':
    import multiprocessing as mp
    from joblib import Parallel, delayed
    from tqdm import tqdm
    import numpy as np
    from math import sqrt
    from multiprocessing import Pool
    import pandas as pd
    from functools import partial

    num_cores = mp.cpu_count()

    ###################################################################################################### worka

    Adata = pd.read_csv(r'C:\Users\Salvo\Desktop\IFISC\data\INDEXES_PBW\INDEX_A_PBW.txt', sep=",", header=None,
                        low_memory=False)

    # total no. of counties
    TOT_counties = len(Adata[1][:]) - 1

    # List of total number of words per county
    tot_word = pd.read_csv('C:\\Users\Salvo\Desktop\IFISC\data\\no_words_per_county.csv', sep=",", low_memory=False,
                           index_col=[0])
    print(tot_word)
    # How many counties you want to analyse? Change first number
    # n_counties = 210 + 1
    # S = 4
    which_county = np.arange(1, 20 + 1)
    # number of words of the i-th county
    word = np.zeros(len(which_county))
    # indx o i-th county
    c = which_county
    print('ciao')

    # Initializing variables
    x = np.arange(1, 70000)
    for i in range(len(which_county)):
        word[i] = int(tot_word.iloc[c[i], 1])
        b = int(word[i])

        y = [[0 for indiceacaso in range(1)] for indiceacaso2 in which_county]
        # cdf_y = [[0 for indiceacaso in range(70000)] for indiceacaso2 in which_county]

    for i in range(0, len(which_county)):
        word[i] = int(tot_word.iloc[c[i], 1])
        print(word[i])
        print(Adata[1][c[i]])
        print(y[i])

        # IMPORTING DATA
        filename = "C:\\Users\Salvo\Desktop\IFISC\data\\all_counties\\row_frequencies\\frequencies_" + str(
            Adata[1][c[i]]) + "_(" + str(c[i]) + ")" + ".txt"
        dat = pd.read_csv(filename, sep=" ", header=None)
        # print(dat)

        b = int(word[i])

        # Real normalized data
        # PMF
        y[i] = dat.iloc[:, 0]
        y[i] = - np.sort(-y[i])
        y[i] = y[i] / np.sum(y[i])
        print('Are the data normalized ?:', np.sum(y[i]))
        print(y[i])
    grid = which_county
    jd_D = np.zeros((len(which_county), len(which_county)))


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

        return js_d


    # c = iterable

    # Define an output queue
    output = mp.Queue()
    # Setup a list of processes that we want to run
    processes = [mp.Process(target=js_distance, args=(y[1], y[j])) for j in range(20)]

    # Run processes
    for p in processes:
        p.start()

    # Exit the completed processes
    for p in processes:
        p.join()


# Get process results from the output queue
results = [output.get() for p in processes]

print(results)


##################################################à##################################################à



'''
if __name__ == '__main__':
    pool = Pool()                         # Create a multiprocessing Pool
    pool.map(js_distance(), data_inputs)
'''

'''
def main():
    iterable = [1, 2, 3, 4, 5]
    pool = multiprocessing.Pool()
    a = "hi"
    b = "there"
    func = partial(js_distance(), a, b)
    pool.map(func, iterable)
    pool.close()
    pool.join()

if __name__ == "__main__":
    main()
'''

'''

def square(x):
    # calculate the square of the value of x
    return x*x

if __name__ == '__main__':

    # Define the dataset
    dataset = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

    # Output the dataset
    print ('Dataset: ' + str(dataset))

    # Run this with a pool of 5 agents having a chunksize of 3 until finished
    agents = 5
    chunksize = 3
    with Pool(processes=agents) as pool:
        result = pool.map(square, dataset, chunksize)

    # Output the result
    print ('Result:  ' + str(result))
'''