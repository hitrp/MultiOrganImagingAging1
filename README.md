# MultiOrganImagingAging
Multi-organ aging clocks characterized by imaging

## 1.Modeling of imaging-based biological age

```R
# modeling of imaging-based biological age for the seven organs with LASSO model.
Rscript age_modeling_optimalPara.R
# the same as above, but stratified by sex.
Rscript age_modeling_optimalPara_stratified.R
```

## 2. Survival analysis

```R
# evaludate the influence of organ age gap on incident risk of diseases and mortality.
Rscript survivalAnalysis.R
# the same as above, but regress out the effect of chronological age, sex and scanning site prior to survival analysis.
Rscript survivalAnalysis_controlSite.R
# the same as above, but group organ-specific age gap into quartiles and take xxx as discrete predictors.
Rscript survivalAnalysis_discrete.R
```

## 3.Association with plasma proteins and blood biomarkers

```R
# Protein
# get the intersection of proteomic and imaging-based age gap data, and remove the effect of chronological age and sex from proteomic data before association analysis.
matlab -nojvm -r processing_metabolism.m
# perform assocation analysis between proteomic and age gap data while controlling for covariates.
Rscript regression.R
# the same as above, but with additional covariate, scanning site.
Rscript regression_withSite.R

# Blood biomarker
# get the intersection of blood biomarker and imaging-based age gap data, and remove the effect of chronological age and sex from blood biomarker data before association analysis.
matlab -nojvm -r processing_proteomic.m
# perform association analysis between blood biomarker and age gap data while controlling for covariates.
Rscript regression_metabolism.R
# the same as above, but with additional covariate, scanning site.
Rscript regression_metabolism_withSite.R
```

## 4.PheWAS

```R
# identify modifiable factors for imaging-based age gap of eye, detailed description of the parameters could be found at https://github.com/MRCIEU/PHESANT.
bash PHESANT_processing_0visit.sh
# identify modifiable factors for imaging-based age gaps of the other six organs, detailed description of the parameters could be found at https://github.com/MRCIEU/PHESANT.
bash PHESANT_processing_2visit.sh
```

## 5. Prediction

```python
# 1. Generate an index for 10-fold cross validation
python merge_index.py
# 2. Training
# The data is trained using lightgbm and divided into three files depending on the situation, namely features only, covariates only and features plus covariates.
# The input parmaters are, form left to right, the location of the source data, the location of the results, and the column's name representing the age(two options, `Age` and `age2`).
python straight_lgbm_feature.py LiangYing_v4_1 trial_v4_1 age2
python straight_lgbm_cov.py LiangYing_v4_1 trial_v4_1 age2
# 3. Evaluation
# From Training, the risk score of each eid can be obtained, and based on this value, AUC and other metrics predicted by the model can be calculated.
# The input parameters are, from left to right, the location of the risk score(a directory under the third parameter), the location of the evaluation results, and the location of the results(same as the second parameter in Training).
python Evalution.py RiskScore_cov Evaluation_cov trial_v4_1
python Evalution.py RiskScore_feature Evaluation_feature trial_v4_1
```

## 6. Plotting

```R
# the code for plotting figures in the main text.
Rscript Plot_code_summary.R
```

