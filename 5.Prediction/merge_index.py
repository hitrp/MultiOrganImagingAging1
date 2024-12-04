import numpy as np
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed
import os
from sklearn.model_selection import train_test_split, StratifiedKFold, KFold

dpath = "/public/share/tmp/LiangY/"
dpath ="/public/share/tmp/LiangYing_v2/"
dpath ="/public/share/tmp/LiangYing_v3/"

items = os.listdir(dpath)
tgt_name_lst = []
for i in items:
    if i.endswith('.csv'):
        name = i[:-4]  # remove suffix of the file
        tgt_name_lst.append(name)

eid_all = pd.DataFrame()
for tgt_name in tqdm(tgt_name_lst):
    eid_df = pd.read_csv(dpath + tgt_name + '.csv')
    result = pd.DataFrame({'eid': eid_df['eid']})
    eid_all = pd.concat([eid_all, result], axis=0).drop_duplicates().reset_index(drop=True)
    
kf = KFold(n_splits=10, shuffle=True, random_state=42)
for i, (_, test_index) in enumerate(kf.split(eid_all)):
    eid_all.loc[eid_all.index[test_index], 'Split'] = i
#result['Split'] = result['Split'].fillna(np.nan).astype('Int64')
eid_all.to_csv("/public/share/tmp/ly_for_rp_swj/trial4/Index.csv",index=False)
