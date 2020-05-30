import numpy as np
import pandas as pd
import time
import myranda as myr
import matplotlib.pyplot as plt


def NameCheck(a, b):
    check = True
    for j in range(a.shape[0]):
        if a[j] not in b:
            print(a[j])
            check = False
    if check:
        print('Yuppie! I dataset coincidono! :)')


pop_tags = ['Spain', 'Parliamentary']
party_tag = 'Spain'
ele_tags = ['Spain', 'parliament']

nation = 'spa'
realistic = False
pseudo_pop = 10000


ele_ds = pd.read_csv('./datasets/view_election.csv', header=0, encoding='iso-8859-1')
party_ds = pd.read_csv('./datasets/view_party.csv', header=0, encoding='iso-8859-1')
pop_ds = pd.read_csv('./datasets/population_dataset.csv', header=0, encoding='utf-8')

pop_ds = pop_ds[pop_ds['Country'] == pop_tags[0]]
pop_ds = pop_ds[pop_ds['Election type'] == pop_tags[1]]
total_pop = pop_ds['Registration'].values.astype(int)
voting_pop = pop_ds['Total vote'].values.astype(int)
year_pop = pop_ds['Year'].values.astype(int)
invalid_pop = pop_ds['Invalid votes'].values
for i in range(invalid_pop.shape[0]):
    invalid_pop[i] = float(invalid_pop[i][:-2])*voting_pop[i]//100
invalid_pop = invalid_pop.astype(int)

party_ds = party_ds[party_ds['country_name'] == party_tag]

ele_ds = ele_ds[ele_ds['country_name'] == ele_tags[0]]
ele_ds = ele_ds[ele_ds['election_type'] == ele_tags[1]]
year_set = np.sort(np.array(list(set([ele_ds['election_date'].values[i][:4]
                                 for i in range(ele_ds['election_date'].shape[0])]))).astype(int))
print(party_ds.head(), party_ds.tail())
print(ele_ds.head(), ele_ds.tail())
other_perc = np.zeros((2, year_set.shape[0]))
other_perc[0] = year_set
other_list = []
for i, c in enumerate(ele_ds['party_name_english']):  # percentages!!!
    if np.nan_to_num(ele_ds['left_right'].values[i]) == 0:
        other_perc[1][np.where(other_perc[0] == int(ele_ds['election_date'].values[i][:4]))]\
            += ele_ds['vote_share'].values[i]
        if c not in other_list:
            other_list.append(c)
for i in other_list:
    ele_ds = ele_ds[ele_ds['party_name_english'] != i]
    party_ds = party_ds[party_ds['party_name_english'] != i]
party_names = party_ds['party_name_english'].values
np.append(party_names, 'others')
np.append(party_names, 'invalid')
np.append(party_names, 'non voters')
party_lr = party_ds['left_right'].values.astype(float)
party_space = np.arange(0, party_lr.shape[0] + 3)
NameCheck(ele_ds['party_name_english'].values, party_names)


# mega_party = (year, array_parties_2d) in cui ogni party è caratterizzato da un int
# array_parties_2d = (idx, n_votes, lr_idx)
mega_party = np.zeros((year_set.shape[0], 3, party_space.shape[0]))
for n, y in enumerate(year_set):
    mega_party[n][0] = party_space
    mega_party[n][2][:-3] = party_lr
    mega_party[n][1][-3] = other_perc[1][n]
    for i in range(ele_ds['election_date'].values.shape[0]):
        if int(ele_ds['election_date'].values[i][:4]) == y:
            mega_party[n][1][np.where(party_names == ele_ds['party_name_english'].values[i])]\
                += ele_ds['vote_share'].values[i]

# year_set è crescente, year_pop decrescente!!!
if realistic:
    difficulty = 'realistic'
    for n in range(year_set.shape[0]):
        mega_party[n][1][-2] = invalid_pop[-n-1]
        mega_party[n][1][-1] = total_pop[-n-1] - voting_pop[-n-1]
        for i in range(mega_party[n][1].shape[0] - 2):
            mega_party[n][1][i] = mega_party[n][1][i] * voting_pop[-n-1] // 100
    np.savez('mp_pn_ys_' + nation + '_' + difficulty + '.npz', mega_party, party_names, year_set)
    print('\nFiles salvati con successo :)')
else:
    difficulty = 'easy'
    for n in range(year_set.shape[0]):
        for i in range(mega_party[n][1].shape[0] - 2):
            mega_party[n][1][i] = mega_party[n][1][i] * pseudo_pop // 100
        mega_party[n][1][-3] = 0
        mega_party[n][1][-2] = 0
        mega_party[n][1][-1] = 0
    np.savez('mp_pn_ys_' + nation + '_' + difficulty + '_' + str(pseudo_pop) +
             '.npz', mega_party, party_names, year_set)
    print('\nFiles salvati con successo :)')


# filez = np.load('mp_pn_ys_' + nation + '_' + difficulty + '_' + str(pseudo_pop) + '.npz')
# print(filez['arr_1'])
