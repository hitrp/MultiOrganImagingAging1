import numpy as np
import pandas as pd
import re
import os
import glob
from tqdm import tqdm
from joblib import Parallel, delayed
from lightgbm import LGBMClassifier
from sklearn.calibration import CalibratedClassifierCV
import sys

nb_cpus = 20
my_seed = 2024
fold_id_lst = [i for i in range(10)]

source = sys.argv[1]
dpath =f"/data/{source}/"
trial = sys.argv[2]
Age = sys.argv[3]

my_params = {'n_estimators': 100,
             'max_depth': 15,
             'num_leaves': 10,
             'subsample': 0.7,
             'learning_rate': 0.01,
             'colsample_bytree': 0.7}

def sort_nicely(l):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c.replace("_","")) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )
    return l

def get_pred_testdata(traindf, testdf, omics_f_lst, my_params, my_seed):
    X_train, y_train = traindf[omics_f_lst], traindf.target_y
    X_test, y_test = testdf[omics_f_lst], testdf.target_y
    my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, seed=my_seed)
    my_lgb.set_params(**my_params)
    calibrate = CalibratedClassifierCV(my_lgb, method='isotonic', cv=5)
    calibrate.fit(X_train, y_train)
    testdf['risk_score'] = calibrate.predict_proba(X_test)[:, 1].tolist()
    pred_df = testdf[['eid', 'Split', 'target_y', 'BL2Target_yrs', 'risk_score']]
    return pred_df

def pred_traindata_iter(traindf, omics_f_lst, my_params, my_seed, fold_id):
    train_idx = traindf['Split'].index[traindf['Split'] != fold_id]
    test_idx = traindf['Split'].index[traindf['Split'] == fold_id]
    X_train, y_train = traindf.iloc[train_idx][omics_f_lst], traindf.iloc[train_idx].target_y
    X_test, y_test = traindf.iloc[test_idx][omics_f_lst], traindf.iloc[test_idx].target_y
    my_lgb = LGBMClassifier(objective='binary', metric='auc', is_unbalance=True, seed=my_seed)
    my_lgb.set_params(**my_params)
    calibrate = CalibratedClassifierCV(my_lgb, method='isotonic', cv=5)
    calibrate.fit(X_train, y_train)
    y_pred_lst = calibrate.predict_proba(X_test)[:, 1].tolist()
    eid_lst = traindf.iloc[test_idx].eid.tolist()
    pred_fold_df = pd.DataFrame({'eid': eid_lst, 'risk_score': y_pred_lst})
    return pred_fold_df

tgt_dir_lst = sort_nicely(glob.glob(f'/data/{source}/*.csv'))
tgt_name_lst = [os.path.basename(tgt_dir)[:-4] for tgt_dir in tgt_dir_lst]

bad_tgt_lst = []
for tgt_name in tqdm(tgt_name_lst):
    mydf = pd.read_csv(dpath + tgt_name + '.csv')
    try:
        omics_f_lst = [Age,"Sex"]
        eid_df = pd.read_csv('/data/trial4/Index.csv')
        traindf = pd.merge(eid_df, mydf[['eid', 'target_y']+omics_f_lst], how='inner', on=['eid'])
        pred_fold_df_lst = Parallel(n_jobs=nb_cpus)(delayed(pred_traindata_iter)(traindf, omics_f_lst, my_params, my_seed, fold_id) for fold_id in fold_id_lst)
        y_pred_train = pd.DataFrame()
        for pred_fold_df in pred_fold_df_lst:
            y_pred_train = pd.concat([y_pred_train, pred_fold_df], axis=0)
        y_pred_train = pd.merge(traindf[['eid', 'Split', 'target_y']], y_pred_train, how='inner', on=['eid'])
        y_pred_train.sort_values(by='eid', ascending=True, inplace=True)
        y_pred_train.to_csv(f'/data/{trial}/RiskScore_cov/' + tgt_name + '.csv',index=False)
    except:
        bad_tgt_lst.append(tgt_name)

bad_tgt_df = pd.DataFrame(bad_tgt_lst)
bad_tgt_df.to_csv(f'/data/{trial}/Bad_RS_c_Target.csv', index=False)