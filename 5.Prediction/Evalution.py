import glob
import re
import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import warnings
import random
from sklearn.metrics import roc_auc_score, confusion_matrix
from sklearn.metrics import roc_curve
from joblib import Parallel, delayed
import sys

#my_omics = 'Proteomic'
#group = 'Incident'
nb_cpus = 20
#nb_iters = 1000
my_seed = 2024
fold_id_lst = [i for i in range(10)]
inpath = sys.argv[1]
outpath = sys.argv[2]
trial = sys.argv[3]

#dpath = '/home1/jiayou/Documents/Projects/BloodOmicsPred/'
#dpath = '/Volumes/JasonWork/Projects/BloodOmicsPred/'
dpath = '/public/share/tmp/LiangY/'
dpath ="/public/share/tmp/LiangYing_v2/"
#dpath ="/public/share/tmp/LiangYing_v3/"
#dpath = "/home1/liangying/LightGBM/ServerCode/for_rp/"


def sort_nicely(l):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c.replace("_","")) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )
    return l


def threshold(array, cutoff):
    array1 = array.copy()
    array1[array1 < cutoff] = 0
    array1[array1 >= cutoff] = 1
    return array1


def Find_Optimal_Cutoff(target, predicted):
    fpr, tpr, threshold = roc_curve(target, predicted)
    i = np.arange(len(tpr))
    roc = pd.DataFrame({'tf': pd.Series(tpr - (1 - fpr), index=i), 'threshold': pd.Series(threshold, index=i)})
    roc_t = roc.iloc[(roc.tf - 0).abs().argsort()[:1]]
    return list(roc_t['threshold'])


def get_eval(y_test, pred_prob, cutoff):
    pred_binary = threshold(pred_prob, cutoff)
    tn, fp, fn, tp = confusion_matrix(y_test, pred_binary).ravel()
    acc = (tp + tn) / (tp + tn + fp + fn)
    sens = tp / (tp + fn)
    spec = tn / (tn + fp)
    prec = tp / (tp + fp)
    Youden = sens + spec - 1
    f1 = 2 * prec * sens / (prec + sens)
    auc = roc_auc_score(y_test, pred_prob)
    evaluations = np.round((auc, acc, sens, spec, prec, Youden, f1), 5)
    evaluations = pd.DataFrame(evaluations).T
    evaluations.columns = ['AUC', 'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'Youden-index', 'F1-score']
    return evaluations
    
'''
def convert_output(result_df):
    result_df = result_df.T
    result_df['Median'] = result_df.median(axis=1)
    result_df['LBD'] = result_df.quantile(0.025, axis=1)
    result_df['UBD'] = result_df.quantile(0.975, axis=1)
    output_lst = []
    for i in range(7):
        output_lst.append('{:.3f}'.format(result_df['Median'][i]) + ' [' +
                          '{:.3f}'.format(result_df['LBD'][i]) + ' - ' +
                          '{:.3f}'.format(result_df['UBD'][i]) + ']')
    result_df['output'] = output_lst
    series = result_df['output']
    return series.to_frame().T

def get_iter_output(mydf, gt_col, y_pred_col, opt_ct, my_iter):
    tmp_random = np.random.RandomState(my_iter)
    bt_idx = tmp_random.choice(range(len(mydf)), size=len(mydf), replace=True)
    mydf_bt = mydf.copy()
    mydf_bt = mydf_bt.iloc[bt_idx, :]
    mydf_bt.reset_index(inplace=True, drop=True)
    y_test_bt = mydf_bt[gt_col]
    eval_iter = get_eval(y_test_bt, mydf_bt[y_pred_col], opt_ct)
    return eval_iter
'''

def get_iter_output(mydf, gt_col, y_pred_col, opt_ct):
    # 直接使用原始数据，不进行抽样
    y_test_bt = mydf[gt_col]
    # 评估模型
    eval_iter = get_eval(y_test_bt, mydf[y_pred_col], opt_ct)
    return eval_iter
    
tgt_dir_lst = sort_nicely(glob.glob(f'/public/share/tmp/ly_for_rp_swj/{trial}/RiskScore_cov/*.csv'))
tgt_name_lst = [os.path.basename(tgt_dir)[:-4] for tgt_dir in tgt_dir_lst]
#tgt_dir_lst = pd.read_csv("/public/share/tmp/ly_for_rp_swj/trail1/Bad_FS_Target.csv")
#tgt_name_lst = tgt_dir_lst["0"].tolist()
'''
tgt_name_lst=[
"liver_M3_0.2_diseasePreidctionData",
"liver_M4_0.2_diseasePreidctionData",
"liver_N0_0.2_diseasePreidctionData",
"liver_N1_0.2_diseasePreidctionData",
"liver_X0_0.2_diseasePreidctionData",
"liver_X1_0.2_diseasePreidctionData",
"liver_X2_0.2_diseasePreidctionData",
"liver_X3_0.2_diseasePreidctionData",
"liver_X4_0.2_diseasePreidctionData"
]
'''
for tgt_name in tqdm(tgt_name_lst):
    tgt_pred_df = pd.read_csv(f'/public/share/tmp/ly_for_rp_swj/{trial}/'+inpath+'/'+tgt_name+'.csv')
    all_eval_lst = []
    for fold_id in fold_id_lst:
        idx = tgt_pred_df['Split'].index[tgt_pred_df['Split'] == fold_id]
        pred_df = tgt_pred_df.iloc[idx]
        unique_classes = pred_df['target_y'].nunique()
        if unique_classes < 2:
            print(f"Warning: Only one class present in fold {fold_id}. Skipping AUC calculation.")
            continue  # 跳过该 fold，避免计算 AUC 错误
        opt_ct = Find_Optimal_Cutoff(pred_df.target_y, pred_df.risk_score)[0]
        iter_eval_lst = get_iter_output(pred_df, 'target_y', 'risk_score', opt_ct)
        all_eval_lst.append(iter_eval_lst)
    all_eval_lst = np.squeeze(all_eval_lst)
    res_eval_df = pd.DataFrame(all_eval_lst, columns=['AUC', 'Accuracy', 'Sensitivity', 'Specificity', 'Precision', 'Youden-index', 'F1-score'])
    mean_res_eval = res_eval_df.mean().round(4)
    mean_res_eval_df = pd.DataFrame(mean_res_eval).T
    mean_res_eval_df.index = ['Mean']
    '''
    eval_df = pd.DataFrame()
    for iter_eval in iter_eval_lst:
        eval_df = pd.concat([eval_df, iter_eval], axis=0)
        print(eval_df)
    res_eval = convert_output(eval_df)
    '''
    mean_res_eval_df.to_csv(f'/public/share/tmp/ly_for_rp_swj/{trial}/'+outpath+'/'+tgt_name+'.csv', index=True)
