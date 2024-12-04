# module R and matlab first
getsampleIntersection <- function(x,y,c,field,excludePatients=FALSE){
  if(!excludePatients){
    # get the intersection of the three matrix by certain column
    temp <- intersect(x[,field],y[,field])
    intersection_field <- intersect(temp, c[,field])
    rownames(x) <- x[,field]
    rownames(y) <- y[,field]
    rownames(c) <- c[,field]
    x[,field] <- NULL
    intersection_x <- as.matrix(x[intersection_field,])
    y[,field] <- NULL
    intersection_y <- as.matrix(y[intersection_field,])
    c[,field] <- NULL
    intersection_c <- as.matrix(c[intersection_field,])

  }else{
    # get the intersection of the three matrix by certain column
    temp <- intersect(x[,field],y[,field])
    intersection_field <- intersect(temp, c[,field])
    patients <- read.table('/data/excludingSubjs.csv',sep=',')
    intersection_field <- setdiff(intersection_field, patients[,1])
    rownames(x) <- x[,field]
    rownames(y) <- y[,field]
    rownames(c) <- c[,field]
    x[,field] <- NULL
    intersection_x <- as.matrix(x[intersection_field,])
    y[,field] <- NULL
    intersection_y <- as.matrix(y[intersection_field,])
    c[,field] <- NULL
    intersection_c <- as.matrix(c[intersection_field,])
  }
  return(list(intersection_x,intersection_y,intersection_c))
}

inverseNorm_transformation <- function(data){
  library(RNOmni)
  return(apply(data,2,function(x)RankNorm(x)))
}
organs <- c('Brain_GM','Brain_WM','heart','boneSub','liver','kidney','Pancreas','OCT')
for(i in seq(1,length(organs))){
organ <- organs[i]
# relationship between Meta and age gap
X <- read.table(paste0('/data/MetabolicAsso/imputated_median_scaled_bl_intersection_',organ,'Age_original_INT_AgeSexremoved_metabolism.csv'),sep = ',',header=FALSE)
colnames(X)[colnames(X)=='V1'] <- 'eid'
X_orig <- as.matrix(X[,2:dim(X)[2]])
Y <- read.table(paste0('/data/MetabolicAsso/imputated_median_scaled_bl_metabolism_intersection_',organ,'Age_original.csv'),sep = ',',header=TRUE)
Y_orig <- as.matrix(Y[,2:dim(Y)[2]])
# check if the subjs in same order
X_order <-  read.table(paste0('/data/MetabolicAsso/imputated_median_scaled_bl_metabolism_intersection_',organ,'Age_original_subjs.csv'),sep = ',',header=TRUE)
Y_order <- read.table(paste0('/data/MetabolicAsso/imputated_median_scaled_bl_metabolism_intersection_',organ,'Age_original_subjs.csv'),sep = ',',header=TRUE)
identical(X_order,Y_order)
eid <- X_order[,1]
# load covariates
cov <- read.table('/data/cova_data_more_withPC_all2visits.csv',sep = ',',header=TRUE)
cov <- cov[!is.na(cov$smoking2) & !is.na(cov$drinking2) & !is.na(cov$Ethnicity) & cov$smoking2!=-3 & cov$drinking2!=-3 & cov$Ethnicity!=0,]
cov_orig <- merge(table(eid),cov,by='eid',all.x=TRUE)
# check for nan values
if(organ=='OCT'){
cov_orig <- cov_orig[,c('eid','age0','Sex','eTIV','Site1','Site2','Education','Townsend_index','BMI0','smoking0','drinking0','Ethnicity')]
cov_orig <- cov_orig[!apply(cov_orig,1,function(x){any(is.nan(x))}),]
cov_orig <- cov_orig[!apply(cov_orig,1,function(x){any(is.na(x))}),]
cov_orig <- fastDummies::dummy_cols(cov_orig, select_columns = "drinking0", remove_first_dummy = TRUE)
cov_orig <- fastDummies::dummy_cols(cov_orig, select_columns = "smoking0", remove_first_dummy = TRUE)
cov_orig <- fastDummies::dummy_cols(cov_orig, select_columns = "Ethnicity", remove_first_dummy = TRUE)
drinkingCol <- colnames(cov_orig)[grepl('drinking0_',colnames(cov_orig))]
smokingCol <- colnames(cov_orig)[grepl('smoking0_',colnames(cov_orig))]
ethnicityCol <- colnames(cov_orig)[grepl('Ethnicity_',colnames(cov_orig))]
cov_orig <- cov_orig[,c('eid','age0','Sex','eTIV','Site1','Site2','Education','Townsend_index','BMI0',smokingCol,drinkingCol,ethnicityCol)]
}else{
cov_orig <- cov_orig[,c('eid','age2','interval','Sex','eTIV','Site1','Site2','Education','Townsend_index','BMI2','smoking2','drinking2','Ethnicity')]
cov_orig <- cov_orig[!apply(cov_orig,1,function(x){any(is.nan(x))}),]
cov_orig <- cov_orig[!apply(cov_orig,1,function(x){any(is.na(x))}),]
cov_orig <- fastDummies::dummy_cols(cov_orig, select_columns = "drinking2", remove_first_dummy = TRUE)
cov_orig <- fastDummies::dummy_cols(cov_orig, select_columns = "smoking2", remove_first_dummy = TRUE)
cov_orig <- fastDummies::dummy_cols(cov_orig, select_columns = "Ethnicity", remove_first_dummy = TRUE)
drinkingCol <- colnames(cov_orig)[grepl('drinking2_',colnames(cov_orig))]
smokingCol <- colnames(cov_orig)[grepl('smoking2_',colnames(cov_orig))]
ethnicityCol <- colnames(cov_orig)[grepl('Ethnicity_',colnames(cov_orig))]
cov_orig <- cov_orig[,c('eid','age2','interval','Sex','eTIV','Site1','Site2','Education','Townsend_index','BMI2',smokingCol,drinkingCol,ethnicityCol)]
}
ageGapNames=organ
MetaNames=read.table('/data/MetabolicAsso/imputated_Metabolicdata_bl_metaNames.csv', sep = ',', header=TRUE)
MetaNames=as.character(MetaNames[,1])
list_samples <- getsampleIntersection(X,Y,cov_orig,'eid',FALSE)
X <- list_samples[[1]]
Y <- list_samples[[2]]
Y <- inverseNorm_transformation(Y)
C <- list_samples[[3]]
modeltype <- 'gaussian'
modality <- organ
prefix <- paste0('imputated_median_scaled_bl_intersection_',organ,'_allremoved')
rownames(cov) <- cov$eid
tempout <- cbind(rownames(C),Y)
colnames(tempout) <- c('eid',ageGapNames)
write.table(tempout,paste0('/data/MetabolicAsso/association_final',organ,'Data_metabolism.csv'),row.names = F,sep=',')
demography_subjs <- rownames(C)
if(organ=='OCT'){
  demography_ageGap <- cov[demography_subjs,c('eid','age0','Sex','eTIV','Site1','Site2','Education','Townsend_index','BMI0','smoking0','drinking0','Ethnicity')]
}else{
  demography_ageGap <- cov[demography_subjs,c('eid','age2','interval','Sex','eTIV','Site1','Site2','Education','Townsend_index','BMI2','smoking2','drinking2','Ethnicity')]
}
demography_ageGap$Group <- rep(organ,nrow(demography_ageGap))
saveRDS(demography_ageGap,paste0('/data/MetabolicAsso/demography_',organ,'_metabolism.rds'))


## model linear regression between multiple Metas and age gap
getCategories <- function(inputMetas){
  assay=read.table('/data/olink_assay.dat', sep='\t', header=TRUE)
  list_assay <- list()
  list_assay$oncology <- assay[grepl('Oncology',assay[,3]),1]
  list_assay$neurology <- assay[grepl('Neurology',assay[,3]),1]
  list_assay$cardioMet <-  assay[grepl('Cardiometabolic',assay[,3]),1]
  list_assay$inflammation <- assay[grepl('Inflammation',assay[,3]),1]
  intersection <- sapply(list_assay,function(x) length(intersect(x, inputMetas)))
  return(intersection)
}
# Define the function to perform linear regression
perform_regression <- function(i, j, modeltype='gaussian',modality) {
  y <- Y[, j]
  x <- X[, i]
  if(modeltype=='gaussian'){
    # BMI2 is not included, it's colinear with eTIV and reduce too much asssociations
    if (modality=='OCT'){
      model <- glm(as.formula(paste0('y~x+age0+Sex+Education+Townsend_index+',paste(smokingCol,sep='+'),'+',paste(drinkingCol,sep='+'),'+',paste(ethnicityCol,sep='+'))), cbind(data.frame(x,y),C),family=gaussian)
    }else{
      model <- glm(as.formula(paste0('y~x+age2+interval+Sex+Education+Townsend_index+',paste(smokingCol,sep='+'),'+',paste(drinkingCol,sep='+'),'+',paste(ethnicityCol,sep='+'))), cbind(data.frame(x,y),C),family=gaussian)
    }
  }
  result <- summary(model)
  return(paste(as.character(result$coefficients[2, 1]),as.character(result$coefficients[2,2]),as.character(result$coefficients[2,3]),as.character(result$coefficients[2,4]),sep='_'))
}
library(foreach)
library(doParallel)
# Perform linear regression for each combination of columns using foreach
num_cores <- detectCores()
registerDoParallel(5)  # Adjust the number of cores as needed
results <- foreach(i = 1:ncol(X), .export = c("X", "Y", 'C','modeltype','modality'), .combine = rbind) %:%
  foreach(j = 1:ncol(Y), .export = c("X", "Y",'C', 'modeltype','modality'), .combine = c) %dopar% {
    perform_regression(i, j,modeltype,modality)
  }
# Apply the function to each element in the matrix
results_beta <- as.matrix(apply(results, c(1, 2), function(x) as.numeric(strsplit(x,'_')[[1]][1])))
results_se <- as.matrix(apply(results, c(1, 2), function(x) as.numeric(strsplit(x,'_')[[1]][2])))
results_tvalue <- as.matrix(apply(results, c(1, 2), function(x) as.numeric(strsplit(x,'_')[[1]][3])))
results_pvalue <- as.matrix(apply(results, c(1, 2), function(x) as.numeric(strsplit(x,'_')[[1]][4])))
padj_BH_overall_vector <- p.adjust(as.vector(results_pvalue), method='BH')
padj_BH_overall <- matrix(padj_BH_overall_vector,nrow = nrow(results_pvalue))
padj_BH_oneMeta <- t(apply(results_pvalue, 1, function(x) p.adjust(x, method='BH')))
padj_Bonf_oneMeta <- t(apply(results_pvalue, 1, function(x) p.adjust(x, method='bonferroni')))
padj_BH_oneageGap <- apply(results_pvalue, 2, function(x) p.adjust(x, method='BH'))
padj_Bonf_oneageGap <- apply(results_pvalue, 2, function(x) p.adjust(x, method='bonferroni'))
statistics_result <- function(pmatrix, dim_index, threshold=0.05){
  # get the number of significant relationships for each Meta
  sig_bool <- pmatrix < threshold
  temp <- results_beta
  temp[!sig_bool]=0
  if(dim_index==1){
    count_sig <- rowSums(sig_bool)
    count_sig_pos <- rowSums(temp>0)
    count_sig_neg <- rowSums(temp<0)
    # get sig agegaps for each Meta
    result_sig <- apply(sig_bool,1,function(x){ageGapNames[x]})
    beta_sig <- apply(temp*sig_bool,1,function(x) x[x!=0])
    se_sig <- apply(results_se*sig_bool,1,function(x) x[x!=0])
    pvalue_sig <- apply(results_pvalue*sig_bool,1,function(x) x[x!=0])
    adjpvalue_sig <- apply(pmatrix*sig_bool,1,function(x) x[x!=0])
  }else if(dim_index==2){
    count_sig <- colSums(sig_bool)
    count_sig_pos <- colSums(temp>0)
    count_sig_neg <- colSums(temp<0)
    # get sig Metas for each column
    result_sig <- apply(sig_bool,2,function(x){MetaNames[x]})
    beta_sig <- apply(temp*sig_bool,2,function(x){x[x!=0]})
    se_sig <- apply(results_se*sig_bool,2,function(x) x[x!=0])
    pvalue_sig <- apply(results_pvalue*sig_bool,2,function(x) x[x!=0])
    adjpvalue_sig <- apply(pmatrix*sig_bool,2,function(x) x[x!=0])
  }
  names(count_sig) <- NULL
  return(list(count_sig, count_sig_pos, count_sig_neg, result_sig,beta_sig,se_sig, pvalue_sig, adjpvalue_sig))
}
count_sig_oneMeta_overall_BH <- statistics_result(padj_BH_overall, 1)[[1]]
count_sig_oneMeta_overall_BH_pos <- statistics_result(padj_BH_overall, 1)[[2]]
count_sig_oneMeta_overall_BH_neg <- statistics_result(padj_BH_overall, 1)[[3]]
agegaps_sig_oneMeta_overall_BH <- statistics_result(padj_BH_overall, 1)[[4]]
agegaps_sig_oneMeta_overall_BH_beta <- statistics_result(padj_BH_overall, 1)[[5]]
agegaps_sig_oneMeta_overall_BH_se <- statistics_result(padj_BH_overall, 1)[[6]]
agegaps_sig_oneMeta_overall_BH_pvalue <- statistics_result(padj_BH_overall, 1)[[7]]
agegaps_sig_oneMeta_overall_BH_adjpvalue <- statistics_result(padj_BH_overall, 1)[[8]]
sig_index_oneMeta_overall_BH <- which(count_sig_oneMeta_overall_BH>0)
count_sig_oneageGap_overall_BH <- statistics_result(padj_BH_overall, 2)[[1]]
count_sig_oneageGap_overall_BH_pos <- statistics_result(padj_BH_overall, 2)[[2]]
count_sig_oneageGap_overall_BH_neg <- statistics_result(padj_BH_overall, 2)[[3]]
Metas_sig_oneageGap_overall_BH <- statistics_result(padj_BH_overall, 2)[[4]]
Metas_sig_oneageGap_overall_BH_beta <- statistics_result(padj_BH_overall, 2)[[5]]
Metas_sig_oneageGap_overall_BH_se <- statistics_result(padj_BH_overall, 2)[[6]]
Metas_sig_oneageGap_overall_BH_pvalue <- statistics_result(padj_BH_overall, 2)[[7]]
Metas_sig_oneageGap_overall_BH_adjpvalue <- statistics_result(padj_BH_overall, 2)[[8]]
sig_index_oneageGap_overall_BH <- which(count_sig_oneageGap_overall_BH>0)

# show distribution across categories
sig_oneMeta_overall_BH <- data.frame(count=count_sig_oneMeta_overall_BH[sig_index_oneMeta_overall_BH], MetaName=MetaNames[sig_index_oneMeta_overall_BH])
sig_oneMeta_overall_BH_all <- data.frame(count=count_sig_oneMeta_overall_BH, count_pos=count_sig_oneMeta_overall_BH_pos, count_neg=count_sig_oneMeta_overall_BH_neg, MetaName=MetaNames)
sig_oneMeta_overall_BH <- sig_oneMeta_overall_BH[order(sig_oneMeta_overall_BH$count, decreasing = TRUE), ]
sig_oneMeta_overall_BH_category <- getCategories(sig_oneMeta_overall_BH[,'MetaName'])
sig_oneageGap_overall_BH <- data.frame(count=count_sig_oneageGap_overall_BH[sig_index_oneageGap_overall_BH], agegapName=ageGapNames[sig_index_oneageGap_overall_BH])
sig_oneageGap_overall_BH <- sig_oneageGap_overall_BH[order(sig_oneageGap_overall_BH$count, decreasing = TRUE), ]
sig_oneageGap_overall_BH_all <- data.frame(count=count_sig_oneageGap_overall_BH, count_pos=count_sig_oneageGap_overall_BH_pos, count_neg=count_sig_oneageGap_overall_BH_neg,agegapName=ageGapNames)


my_data_list <- list(sig_oneMeta_overall_BH, sig_oneMeta_overall_BH_all, sig_oneageGap_overall_BH,sig_oneageGap_overall_BH_all, agegaps_sig_oneMeta_overall_BH,agegaps_sig_oneMeta_overall_BH_beta,agegaps_sig_oneMeta_overall_BH_se,agegaps_sig_oneMeta_overall_BH_pvalue,agegaps_sig_oneMeta_overall_BH_adjpvalue,Metas_sig_oneageGap_overall_BH,Metas_sig_oneageGap_overall_BH_beta,Metas_sig_oneageGap_overall_BH_se,Metas_sig_oneageGap_overall_BH_pvalue,Metas_sig_oneageGap_overall_BH_adjpvalue,ageGapNames,MetaNames,sig_oneMeta_overall_BH_category, results_beta, results_se, results_tvalue, results_pvalue, padj_BH_overall)
names(my_data_list) <- c('sig_oneMeta_overall_BH','sig_oneMeta_overall_BH_all','sig_oneageGap_overall_BH','sig_oneageGap_overall_BH_all','agegaps_sig_oneMeta_overall_BH','agegaps_sig_oneMeta_overall_BH_beta','agegaps_sig_oneMeta_overall_BH_se','agegaps_sig_oneMeta_overall_BH_pvalue','agegaps_sig_oneMeta_overall_BH_adjpvalue','Metas_sig_oneageGap_overall_BH','Metas_sig_oneageGap_overall_BH_beta','Metas_sig_oneageGap_overall_BH_se','Metas_sig_oneageGap_overall_BH_pvalue','Metas_sig_oneageGap_overall_BH_adjpvalue','ageGapNames','MetaNames','sig_oneMeta_overall_BH_category','results_beta','results_se','results_tvalue','results_pvalue','padj_BH_overall')
saveRDS(my_data_list, file = paste('/data/MetabolicAsso/',prefix,"_result_metabolism.rds",sep=''))
save.image(file = paste('/data/MetabolicAsso/',prefix,"_environment_metabolism.RData",sep=''))
}


# actually used to get p of all Metas
organs <- c('Brain_GM','Brain_WM','heart','boneSub','liver','kidney','Pancreas')
for(i in seq(1,length(organs))){
  organ <- organs[i]
  result <- readRDS(paste0('/data/MetabolicAsso/imputated_median_scaled_bl_intersection_',organ,'_allremoved_result_metabolism.rds'))
  temp <- data.frame(result$MetaNames, result$results_beta, result$results_pvalue, result$padj_BH_overall)
  colnames(temp) <- c('MetaName', 'Beta','pvalue', 'adjP') 
  write.table(temp, paste0('/data/MetabolicAsso/imputated_median_scaled_bl_intersection_',organ,'_allremoved_result_metabolism_sigResultBeta.csv'),sep=',',row.names = F)
}
