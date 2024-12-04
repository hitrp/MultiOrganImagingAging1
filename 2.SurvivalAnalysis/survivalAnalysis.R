# survial analysis
library(survival)
library(data.table)
organs <- c('meanOverall','Brain_GM','Brain_WM','heart','kidney','liver','OCT','boneSub','Pancreas')
agingfiles <- sapply(organs,function(x) paste0(x,'_age_prediction_all.csv'))
targetfiles <- list.files('/Targets_RAW/',pattern = '.csv')
covafile <- read.table('/data/cova_data_more_withPC_all2visits.csv',sep=',',header=T)
covadata <- covafile[,c('eid','interval','Sex','Townsend_index','smoking0','smoking2','drinking0','drinking2','BMI0','BMI2','Education','Ethnicity')]
diseaseCode <- read.table('/data/Target_code.csv',sep=',',header=T)
diseaseCode <- diseaseCode[,c('Disease_code','Disease')]
colnames(diseaseCode) <- c('target','Disease')
# cox harzard model, with delta as continuous variable, in all subjs, with more covariates
list_agingfile <- c()
list_targetfile <- c()
list_numberAll <- c()
list_numberTarget <- c()
list_pvalues <- c()
list_zvalues <- c()
list_coef <- c()
list_expcoef <- c()
list_secoef <- c()
list_up <- c()
list_down <- c()
for(i in seq(1,length(agingfiles))){
  organ <- organs[i]
  agingfile <- agingfiles[i]
  agingdata <- read.table(paste0('/data/ageModeling/',agingfile),sep=',',header=T)
  pvalues <- c()
  for(j in seq(1,length(targetfiles))){
    targetfile <- targetfiles[j]
    list_agingfile <- c(list_agingfile, agingfile)
    list_targetfile <- c(list_targetfile, targetfile)
    # convert variables
    if(organ!='OCT'){
      covadata_used <- covadata[,c('eid','interval','Sex','Townsend_index','smoking2','drinking2','BMI2','Education','Ethnicity')]
      covadata_used <- covadata_used[!is.na(covadata_used$smoking2) & !is.na(covadata_used$drinking2) & covadata_used$smoking2!=-3 & covadata_used$drinking2!=-3,]
      covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "drinking2", remove_first_dummy = TRUE)
      covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "smoking2", remove_first_dummy = TRUE)
      tempdata <- merge(agingdata, covadata_used, by=c('eid'))
    }else{
      covadata_used <- covadata[,c('eid','interval','Sex','Townsend_index','smoking0','drinking0','BMI0','Education','Ethnicity')]
      covadata_used <- covadata_used[!is.na(covadata_used$smoking0) & !is.na(covadata_used$drinking0) & covadata_used$smoking0!=-3 & covadata_used$drinking0!=-3,]
      covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "drinking0", remove_first_dummy = TRUE)
      covadata_used <- fastDummies::dummy_cols(covadata_used, select_columns = "smoking0", remove_first_dummy = TRUE)
      tempdata <- merge(agingdata, covadata_used, by=c('eid'))
    }
    targetdata <- read.csv(paste0('/Targets_RAW/',targetfile),sep=',',header=T)
    useddata <- merge(tempdata, targetdata, by=c('eid'))
    if(organ=='OCT'){
      useddata$time <- useddata$BL2Target_yrs
    }else{
      useddata$time <- useddata$BL2Target_yrs-useddata$interval 
    }
    useddata <- useddata[useddata$time>0,]
    list_numberAll <- c(list_numberAll, nrow(useddata))
    list_numberTarget <- c(list_numberTarget, nrow(useddata[useddata$target_y==1,]))
    # survival analysis
    if(organ!='OCT'){
      cleaned_data <- useddata[,c('Age','delta_corrected','Sex','target_y','time','smoking2_1','smoking2_2','drinking2_1','drinking2_2','BMI2','Education')]
      colnames(cleaned_data) <- c('Age','delta_corrected','Sex','status','years','smoking2_1','smoking2_2','drinking2_1','drinking2_2','BMI2','Education')
      cleaned_data$delta_corrected <- scale(cleaned_data$delta_corrected)
      cox_fit <- coxph(Surv(years, status) ~ delta_corrected+Sex+Age+smoking2_1+smoking2_2+drinking2_1+drinking2_2+BMI2+Education, data=cleaned_data)
    }else{
      cleaned_data <- useddata[,c('Age','delta_corrected','Sex','target_y','time','smoking0_1','smoking0_2','drinking0_1','drinking0_2','BMI0','Education')]
      colnames(cleaned_data) <- c('Age','delta_corrected','Sex','status','years','smoking0_1','smoking0_2','drinking0_1','drinking0_2','BMI0','Education')
      cleaned_data$delta_corrected <- scale(cleaned_data$delta_corrected)
      cox_fit <- coxph(Surv(years, status) ~ delta_corrected+Sex+Age+smoking0_1+smoking0_2+drinking0_1+drinking0_2+BMI0+Education, data=cleaned_data)
    }
    result <- summary(cox_fit)$coefficients
    CI_range <- summary(cox_fit)$conf.int[1,c(3,4)]
    # Check for the possibility of infinite coefficients
    if (any(is.infinite(coefficients(cox_fit)))) {
      print("Warning: Model coefficients may be infinite.")
    }
    list_pvalues <- c(list_pvalues, result[1,5])
    list_zvalues <- c(list_zvalues, result[1,4])
    list_coef <- c(list_coef, result[1,1])
    list_expcoef <- c(list_expcoef, result[1,2])
    list_secoef <- c(list_secoef, result[1,3])
    list_up <- c(list_up, CI_range[2])
    list_down <- c(list_down, CI_range[1])
    pvalues <- c(pvalues, result[1,5])
  }
}
resultframe <- data.frame(aging=list_agingfile, target=list_targetfile, numberAll=list_numberAll, numberTarget=list_numberTarget, expcoef=list_expcoef, coef=list_coef, down=unname(list_down), up=unname(list_up), secoef=list_secoef, zvalue=list_zvalues, pvalue=list_pvalues, separateFDR=list_FDR_separate)
resultframe$target <- sapply(resultframe$target, function(x) gsub('.csv','',x))
resultframe <- merge(resultframe, diseaseCode, by='target', all.x = T)
resultframe$Padj_FDR_overall <- p.adjust(resultframe$pvalue,method='BH')
write.table(resultframe, '/data/Survival/allSubjs_cox_moreCova.csv',sep=',',row.names = F)