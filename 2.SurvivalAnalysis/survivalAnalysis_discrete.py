import numpy as np
import math
import pandas as pd
from sklearn.preprocessing import StandardScaler
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import matplotlib as mpl
from lifelines import CoxPHFitter
mpl.rcParams['pdf.fonttype'] = 42



organ_all = ['boneSub', 'Brain_GM', 'Brain_WM', 'heart', 'kidney', 'liver', 'Pancreas','mean']
ill_list = ['A0','A1','A2','C0','C1','C2','C3','C4','C5','C6','D0','D1','E0','E1','E2','F0','F1',
        'F2','F3','F4','F5','G0','G1','G2','G3','G4','H0','H1','H2','H3','I0','I1','I2','I3',
        'I4','I5','I6','I7','J0','J1','J2','J3','K0','K1','K2','K3','L0','L1','M0','M1','M2',
        'M3','M4','N0','N1','X0','X1','X2','X3','X4']

# load the physical data
df_ph = pd.read_csv("/data/cova_data_more_withPC_all2visits.csv",
                    usecols=['eid', 'Sex', 'BMI2', 'smoking2', 'drinking2'])
TDI = pd.read_csv("/data/TDI.csv",
                  usecols=['eid', 'TDI_imp'])
df_ph = pd.merge(df_ph, TDI, how='inner', on=['eid'])
df_ph.rename(columns={'TDI_imp': 'TDI'}, inplace=True)

for organ in organ_all:
    print(organ)
    row_num = 0
    # load the data
    df = pd.read_csv(""+organ+"_age_prediction_all.csv",
                     usecols=['eid', 'delta_corrected', 'Age'])

    df['delta_c'] = df['delta_corrected']
    df['tmp'] = df['delta_c']
    mydf = df

    #normalization
    norm_col = mydf.columns[-2:]
    scaler = StandardScaler()
    mydf[norm_col] = scaler.fit_transform(mydf[norm_col])

    #split the data and mark the label
    mydf.sort_values(by="delta_c", inplace=True, ascending=True)
    mydf = mydf.reset_index(drop=True)

    Q = round(mydf.shape[0] / 4)
    Q_1 = Q - 1
    Q_2 = mydf.shape[0] - Q
    Q_3 = Q_2 - 1

    mydf.loc[:Q_1, 'age_gap_c'] = 0  # 1/4-1
    mydf.loc[Q_2:, 'age_gap_c'] = 1  # all-9241
    mydf.loc[Q:Q_3, 'age_gap_c'] = 2

    #load illness data
    mydf_use = mydf
    age0 = pd.read_csv("/data/cova_data_more_withPC_all2visits.csv",
                       usecols=['eid', 'age0'])
    mydf_use = pd.merge(mydf, age0, how='inner', on=['eid'])
    mydf_use['time'] = mydf_use['Age'] - mydf_use['age0']

    HR = pd.read_csv("HR_all_lianxu.csv")

    for ill in ill_list:
        ill_df = pd.read_csv("/data/Targets_RAW/Targets_RAW/"+ill+".csv",
                            usecols=['eid', 'target_y', 'BL2Target_yrs'])
        ill_df = pd.merge(mydf_use, ill_df, how='inner', on=['eid'])
        ill_df['BL2Target_yrs1'] = ill_df['BL2Target_yrs'] - ill_df['time']
        ill_df.drop(index=ill_df['BL2Target_yrs1'].index[ill_df['BL2Target_yrs1'] < 0], inplace=True)

        ill_ph_df = pd.merge(ill_df, df_ph, how='inner', on=['eid'])

        #cox analyse
        my_formula = "Age + C(Sex) + TDI + BMI2 + C(drinking2) + C(smoking2) + C(age_gap_c)"  # 0 reference  + C(Smoking)   delta_c
        #     my_formula = "Age + C(Sex) + C(age_gap_c)"
        cph = CoxPHFitter()
        if row_num == 17:
            row_num += 1
            continue
        cph.fit(ill_ph_df, duration_col='BL2Target_yrs1', event_col='target_y', formula=my_formula)
        print(row_num)
        HR.loc[row_num,'HR'] = cph.summary['exp(coef)'][-2]
        HR.loc[row_num,'LOW'] = cph.summary['exp(coef) lower 95%'][-2]
        HR.loc[row_num,'UP'] = cph.summary['exp(coef) upper 95%'][-2]
        HR.loc[row_num,'P'] = cph.summary['p'][-2]

        HR.loc[row_num,'all_num'] = ill_ph_df.shape[0]
        HR.loc[row_num,'target_num'] = len(ill_ph_df['target_y'].index[ill_ph_df['target_y'] == 1])
        row_num += 1

        plt.figure(figsize=(8, 8))
        kmf = KaplanMeierFitter()
        j = 0
        # label = ['BA=CA','BA<CA','BA>CA']
        label = ['Q1', 'Q4']
        colors = ['blue', 'red']
        for name, grouped_df in ill_ph_df.groupby('age_gap_c'):
            kmf.fit(grouped_df['BL2Target_yrs1'], grouped_df['target_y'], label=label[j])
            kmf.plot_survival_function(color=colors[j])  # 使用自定义颜色
            j += 1
            if j == 2:
                break
        plt.title('KM Survival Curve')
        plt.xlabel('Time (years)')
        plt.ylabel('Survival Probability')
        plt.legend()
        # plt.show()
        plt.savefig( ""+organ+"_"+ill+"_survival_curve.pdf", format='pdf')

    HR.to_csv(""+organ+"_HR.csv", index=False)

print("finish!")

