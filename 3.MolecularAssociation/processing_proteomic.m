preprocessingMethod = 'imputated_median_scaled';
organs={'Brain_WM','Brain_GM','heart','kidney','liver','boneSub','OCT','Pancreas'};
for organ_index=1:length(organs)
organ = organs{organ_index};
% read prepocessed proteomic 
dataArray = readtable(['/data/' preprocessingMethod '_data_bl.csv'],'Delimiter',',');
dataArray.Properties.RowNames = arrayfun(@(x) num2str(x),dataArray.eid,'UniformOutput',false);
proteinNames = readtable(['/data/' preprocessingMethod '_bl_proteinNames.csv'],'Delimiter',',');

% get the intersection of proteomic and neuroimaging data
ageGap_table = readtable(['/data/ageModeling/',organ,'_age_prediction_all.csv'],'Delimiter',',','ReadVariableNames',true);
ageGap_table = ageGap_table(:,{'eid','delta_corrected'});
regionalName = {'delta_corrected'};
ageGap_table.Properties.RowNames = arrayfun(@num2str,ageGap_table.eid,'UniformOutput',false);

% correct for site and eTIV
cov_data = readtable('/data/cov_data_ukb_20230410.csv');
cov_data.eTIV(cellfun(@(x) isequal(x, 'NA'),cov_data.eTIV)) = {'NaN'};
cov_data.Site1(cellfun(@(x) isequal(x, 'NA'),cov_data.Site1)) = {'NaN'};
cov_data.Site2(cellfun(@(x) isequal(x, 'NA'),cov_data.Site2)) = {'NaN'};
cov_data.eTIV = cellfun(@str2double, cov_data.eTIV);
cov_data.Site1 = cellfun(@str2double, cov_data.Site1);
cov_data.Site2 = cellfun(@str2double, cov_data.Site2);
% filtering for non missing eTIV information
cov_data = cov_data(~any(isnan(cov_data.eTIV),2),:);
cov_data = cov_data(:,{'eid','Age','Sex'});
cov_data.Properties.RowNames = arrayfun(@(x) num2str(x),cov_data.eid,'UniformOutput',false);

% perfrom intersection
intersection_table = innerjoin(dataArray,ageGap_table, "Keys",{'eid'});
intersection_table = innerjoin(intersection_table,cov_data, "Keys",{'eid'});
intersection_subjs = arrayfun(@num2str,intersection_table.eid,'UniformOutput',false);
intersection_ageGap = ageGap_table(intersection_subjs,:);
intersection_protemic = dataArray(intersection_subjs,:);
% remove age and sex effect from proteomic data
% perfrom intersection and remove the effect of age and sex
used_subjs = intersection_subjs;
writetable(table(intersection_subjs),['/data/ProteomicAsso/' preprocessingMethod '_bl_intersection_',organ,'Age_original_subjs.csv'],'Delimiter',',');
writetable(intersection_ageGap,['/data/ProteomicAsso/' preprocessingMethod '_bl_intersection_',organ,'Age_original.csv'],'Delimiter',',');

used_proteomic = intersection_protemic(used_subjs,:);
used_cov = cov_data(used_subjs,:);
C = used_cov{:,2:size(used_cov,2)};
X_before = used_proteomic{:,2:size(used_proteomic,2)};
% impute missing values before regression
encoding = 'UTF-8';
% missing data imputation
% Find the indices of missing values
missingValues = isnan(X_before);
% Replace missing values with the respective column median
for col = 1:size(X_before, 2)
    colMeans = median(X_before(:, col),"omitnan");
    X_before(missingValues(:, col), col) = colMeans;
end

writetable(array2table(X_before),'/temp.csv','Delimiter',',','WriteVariableNames',false);
% perfrom inverse normal transformation with inverseNorm_transformation.R
% inverse normal transformation
% Define the Rscript path and your R script
RscriptPath = 'Rscript';  % Path to Rscript executable
scriptPath = '/script/inverseNorm_transformation.R';  % Your R script
cmd = sprintf('%s %s', RscriptPath, scriptPath);
[status, result]=system(cmd);

X_before = readtable('/temp_INTransform.csv','Delimiter',',');
X_before = table2array(X_before);
C = cellfun(@str2double, C);
% Perform linear regression to estimate the coefficients
B = C\X_before;  % Coefficients matrix (size: number of covariates x number of columns in X)
% Predict the covariate effect
C_effect = C * B;
% Remove the effect of covariates from each column of X
X_corrected = X_before - C_effect;
eid = cellfun(@str2num, intersection_subjs);
writetable(array2table([eid X_corrected]),['/data/ProteomicAsso/' preprocessingMethod '_bl_intersection_',organ,'Age_original_INT_AgeSexremoved_protemic.csv'],'Delimiter',',','WriteVariableNames',false);
end