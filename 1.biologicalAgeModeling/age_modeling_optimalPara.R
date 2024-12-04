library(dplyr)
library(readr)
library(caret)
library(glmnet)
library(Metrics)

# Define paths
data_path <- "/data/data/"
cova_data_path <- "/data/cova_data_more_withPC_all2visits_withheadmotion.csv"
cova_data <- read.csv(cova_data_path)
eid_data_path <- "/data/BA_health_eid.csv"
allSiteInfo <- read.table('/data/AllSiteInfo.csv', sep=',' ,header=T, check.names=F)
HeightInfo <- read.table("/data/AllHeightInfo.csv", sep=',', header=T, check.names=F)
HeightInfo <- HeightInfo[,c('eid','50-0.0')]
colnames(HeightInfo) <- c('eid', 'Height')
# Define organs and their respective data files
organs <- c('Brain_WM', 'Brain_GM', 'heart', 'kidney', 'liver', 'Pancreas', 'OCT', 'boneSub')
organs_data <- c('Brain_WM_clean.csv', 'Brain_GM_clean.csv','heart_clean.csv', 'kidney_clean.csv', 'liver_clean.csv','Pancreas_clean.csv', 'OCT_clean.csv','boneSub_clean.csv')
sink('/data/ageModeling/imagingvisit_ageModeling_log.txt')
########################## model traning and testing on normal controls
for (i in 1:length(organs)) {
  print(organs[i])
  organ <- organs[i]
  organ_data <- organs_data[i]
  
  # Load data
  df <- read.csv(paste0(data_path, organ_data))
  
  # using only normal controls
  eid <- read.csv(eid_data_path)
  df <- inner_join(df, eid['eid'], by = 'eid')
  
  # Remove subjects without complete data
  cnt <- rowSums(!is.na(df))
  df <- df[cnt == ncol(df),]
  col <- colnames(df)
  
  if (organ != 'OCT') {
    age <- read.csv(cova_data_path)
    age <- age[,c('eid','age2')]
    df <- inner_join(df, age, by = 'eid')
    df <- rename(df, Age = age2)
  } else {
    age <- read.csv(cova_data_path)
    age <- age[,c('eid','Age')]
    df <- inner_join(df, age, by = 'eid')
  }
  # remove subjects without real age info
  df <- df[!is.na(df$Age),]
  # get confounding info
  if(organ %in% c('brain','Brain_WM','Brain_GM')){
      df <- merge(df, allSiteInfo[,c('eid','54-2.0')])
      df <- df[df$'54-2.0'!='11028',] # remove the site with too little subjects
      colnames(df)[which(colnames(df) == "54-2.0")] <- "Site"
      df$Site <- as.character(df$Site)
      df <- fastDummies::dummy_cols(df, select_columns = "Site", remove_first_dummy = TRUE)
      SiteCols <- unlist(sapply(colnames(df),function(x) if(grepl('Site_',x)) x),use.names = F)
      TIV <- cova_data
      TIV <- TIV[,c('eid','eTIV','Sex')]
      df <- merge(df, TIV, by = 'eid')
      df <- na.omit(df)
    }else if(organ=='OCT'){
      df <- merge(df, allSiteInfo[,c('eid','54-0.0')])
      df <- df[df$'54-0.0'!='11022',]
      colnames(df)[which(colnames(df) == "54-0.0")] <- "Site"
      df$Site <- as.character(df$Site)
      df <- fastDummies::dummy_cols(df, select_columns = "Site", remove_first_dummy = TRUE)
      SiteCols <- unlist(sapply(colnames(df),function(x) if(grepl('Site_',x)) x),use.names = F)
      df <- merge(df, HeightInfo, by = 'eid')
      Age <- cova_data
      Age <- Age[,c('eid','Sex')]
      df <- merge(df, Age, by = 'eid')
      df <- na.omit(df)
    }else{
      df <- merge(df, allSiteInfo[,c('eid','54-2.0')])
      colnames(df)[which(colnames(df) == "54-2.0")] <- "Site"
      df$Site <- as.character(df$Site)
      df <- fastDummies::dummy_cols(df, select_columns = "Site", remove_first_dummy = TRUE)
      SiteCols <- unlist(sapply(colnames(df),function(x) if(grepl('Site_',x)) x),use.names = F)
      df <- merge(df, HeightInfo, by = 'eid')
      Age <- cova_data
      Age <- Age[,c('eid','Sex')]
      df <- merge(df, Age, by = 'eid')
      df <- na.omit(df)
    }
  table(df$Site)
  # Split data into training and testing sets
  set.seed(12345)
  folds <- createFolds(df$Age, k = 10)
  
  list_mae <- list()
  list_r <- list()
  list_models <- list()
  list_correctionModel <- list()
  list_lambda <- list()
  list_scaler <- list()
  list_plots <- list()
  list_coef <- list()
  list_covModles <- list()
  
  
  alldata <- df
  for (k in 1:10) {
    print(k)
    test_indices <- folds[[k]]
    train_indices <- setdiff(1:nrow(alldata), test_indices)

    x_train_data <- alldata[train_indices,]
    x_test_data <- alldata[test_indices,]
    data_matrix <- x_train_data[,unlist(sapply(colnames(x_train_data),function(x) grepl('^X',x)), use.names = F)]
    data_matrix_test <- x_test_data[,unlist(sapply(colnames(x_test_data),function(x) grepl('^X',x)), use.names = F)]

    residuals_matrix <- matrix(NA, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
    residuals_matrix_test <- matrix(NA, nrow = nrow(data_matrix_test), ncol = ncol(data_matrix_test))
    colnames(residuals_matrix) <- colnames(data_matrix)
    colnames(residuals_matrix_test) <- colnames(data_matrix_test)
    covModels <- list()
    for (cc in 1:ncol(data_matrix)) {
      if(organ %in% c('brain','Brain_WM','Brain_GM')){
        temp_data <- cbind(dependent_variable = data_matrix[, cc], x_train_data[,c('eTIV',SiteCols)])
        testing_data <- cbind(dependent_variable = data_matrix_test[, cc], x_test_data[,c('eTIV',SiteCols)])
        model <- lm(as.formula(paste0('dependent_variable ~ eTIV +',paste(SiteCols,collapse ='+'))), data = temp_data)
      }else{
        temp_data <- cbind(dependent_variable = data_matrix[, cc], x_train_data[,c('Height',SiteCols)])
        testing_data <- cbind(dependent_variable = data_matrix_test[, cc], x_test_data[,c('Height',SiteCols)])
        model <- lm(as.formula(paste0('dependent_variable ~ Height+',paste(SiteCols,collapse = '+'))), data = temp_data)
      }
      residuals_matrix[, cc] <- residuals(model)
      # Predict the effect of covariates on the testing data using the training model
      predicted_values_test <- predict(model, newdata = testing_data)
      # Calculate residuals for the testing data
      residuals_matrix_test[,cc] <- testing_data$dependent_variable - predicted_values_test
      covModels[[cc]] <- model
    }
    list_covModles[[k]] <- covModels
    X_train_fold <- residuals_matrix
    y_train_fold <- alldata[train_indices, 'Age']
    X_test_fold <- residuals_matrix_test
    y_test_fold <- alldata[test_indices, 'Age']
    
    # normalizaion and standardziation
    library(caret)
    normParam <- preProcess(X_train_fold)
    X_train_fold <- predict(normParam, X_train_fold)
    X_test_fold <- predict(normParam, X_test_fold)
    list_scaler[[k]] <- normParam
    
    require(doMC)
    registerDoMC(cores = 10)
    set.seed(1010)
    lasso_reg = glmnet(as.matrix(X_train_fold), y_train_fold, alpha = 1)
    cv_lasso <- cv.glmnet(as.matrix(X_train_fold), y_train_fold, lambda = 10^seq(5, -8, by = -.1), alpha = 1, type.measure = 'mae',nfolds = 10,parallel = T)
    plot(cv_lasso)
    list_plots[[k]] <- recordPlot()
    optimal_lambda <- cv_lasso$lambda.min
    lasso_model = predict(lasso_reg, type = "coefficients", s = optimal_lambda)
    lasso_model_nonzero <- lasso_model[(lasso_model[,1]!= 0),]
    list_lambda[[k]] <- optimal_lambda
    lasso_predict <- predict(lasso_reg, newx = as.matrix(X_test_fold), s = optimal_lambda)
    
    df[test_indices, 'pre_age'] <- as.numeric(lasso_predict)
    
    # Age bias correction
    train_predicted_y <- predict(lasso_reg, newx = as.matrix(X_train_fold), s = optimal_lambda)
    model <- lm(train_predicted_y ~ y_train_fold)
    test <- data.frame(TrueAge = as.numeric(y_test_fold), age_LASSO = as.numeric(lasso_predict))
    test$age_correct_LASSO <- (test$age_LASSO - coef(model)[1]) / coef(model)[2]
    
    df[test_indices, 'delta_corrected'] <- unlist(test$age_correct_LASSO - test$TrueAge)
    df[test_indices, 'pre_age_corrected'] <- test$age_correct_LASSO
    
    list_mae[[k]] <- mae(test$TrueAge, test$age_LASSO)
    list_r[[k]] <- cor(test$TrueAge, test$age_LASSO)
    list_models[[k]] <- lasso_reg
    list_coef[[k]] <- lasso_model
    list_correctionModel[[k]] <- model
  }
  
  best_model_index <- which.min(list_mae)
  best_model <- list_models[[best_model_index]]
  best_correction_model <- list_correctionModel[[best_model_index]]
  best_scaler <- list_scaler[[best_model_index]]
  best_plot <- list_plots[[best_model_index]]
  best_lambda <- list_lambda[[best_model_index]]
  best_coef <- list_coef[[best_model_index]]
  print(list_lambda[[best_model_index]])
  best_covModel <- list_covModles[[best_model_index]]
  
  df <- df[, !(names(df) %in% c(unlist(sapply(colnames(df),function(x) {if(grepl('Site',x)) x}), use.names = F),'eTIV','Height','Sex'))]
  df_HC <- data.frame(df)
  write_csv(df_HC, paste0('/data/ageModeling/', organ, '_age_prediction.csv'))
  save.image(file = paste('/data/ageModeling/',organ,"_AgeModeling_environment.RData",sep=''))

  ##################### Model evaluation
  r <- cor(df$Age, df$pre_age)
  r_corrected <- cor(df$Age, df$pre_age_corrected)
  mae <- mae(df$Age, df$pre_age)
  mae_corrected <- mae(df$Age, df$pre_age_corrected)
  corrdata <- data.frame(df$Age, df$pre_age)
  write.table(corrdata, paste0('/data/ageModeling/', organ, '_age_prediction_evaluation.csv'), sep=',', row.names=F)

  print(r)
  print(r_corrected)
  print(mae)
  print(mae_corrected)

  # Get weights of the predictors
  brain_traits <- data.frame(rownames(best_coef),as.matrix(best_coef))
  colnames(brain_traits) <- c('field','coef')
  write_csv(brain_traits, paste0('/data/ageModeling/', organ, 'Age_LASSO_weights.csv'))


  ##################### Apply the model trained on normal subjects to other subjects (excluding normal controls)
  df <- read.csv(paste0(data_path, organ_data))
  # exclude only normal controls
  eid <- read.csv(eid_data_path)
  df <- df[!(df$eid %in% eid$eid),]
  cnt <- rowSums(!is.na(df))
  df <- df[cnt == ncol(df), ]

  if (organ != 'OCT') {
    age <- read.csv(cova_data_path)
    age <- age[,c('eid','age2')]
    df <- inner_join(df, age, by = 'eid')
    df <- rename(df, Age = age2)
  } else {
    age <- read.csv(cova_data_path)
    age <- age[,c('eid','Age')]
    df <- inner_join(df, age, by = 'eid')
  }
  # remove subjs without real age info
  df <- df[!is.na(df$Age),]
  
  # get confounding info
  if(organ %in% c('brain','Brain_WM','Brain_GM')){
    df <- merge(df, allSiteInfo[,c('eid','54-2.0')])
    df <- df[df$'54-2.0'!='11028',] # remove the site with too little subjects
    colnames(df)[which(colnames(df) == "54-2.0")] <- "Site"
    df$Site <- as.character(df$Site)
    df <- fastDummies::dummy_cols(df, select_columns = "Site", remove_first_dummy = TRUE)
    SiteCols <- unlist(sapply(colnames(df),function(x) if(grepl('Site_',x)) x),use.names = F)
    TIV <- read.csv(cova_data_path)
    TIV <- TIV[,c('eid','eTIV')]
    df <- merge(df, TIV, by = 'eid')
    df <- na.omit(df)
  }else if(organ=='OCT'){
    df <- merge(df, allSiteInfo[,c('eid','54-0.0')])
    df <- df[df$'54-0.0'!='11022',]
    colnames(df)[which(colnames(df) == "54-0.0")] <- "Site"
    df$Site <- as.character(df$Site)
    df <- fastDummies::dummy_cols(df, select_columns = "Site", remove_first_dummy = TRUE)
    SiteCols <- unlist(sapply(colnames(df),function(x) if(grepl('Site_',x)) x),use.names = F)
    df <- merge(df, HeightInfo, by = 'eid')
    df <- na.omit(df)
  }else{
    df <- merge(df, allSiteInfo[,c('eid','54-2.0')])
    colnames(df)[which(colnames(df) == "54-2.0")] <- "Site"
    df$Site <- as.character(df$Site)
    df <- fastDummies::dummy_cols(df, select_columns = "Site", remove_first_dummy = TRUE)
    SiteCols <- unlist(sapply(colnames(df),function(x) if(grepl('Site_',x)) x), use.names = F)
    df <- merge(df, HeightInfo, by = 'eid')
    df <- na.omit(df)
  }
  
  # apply the covariate regression model to other subjs
  x_test_data <- df
  data_matrix_test <- x_test_data[,unlist(sapply(colnames(x_test_data),function(x) grepl('^X',x)), use.names = F)]
  # Initialize a matrix to store residuals
  residuals_matrix_test <- matrix(NA, nrow = nrow(data_matrix_test), ncol = ncol(data_matrix_test))
  colnames(residuals_matrix_test) <- colnames(data_matrix_test)
  for (cc in 1:ncol(data_matrix_test)) {
    if(organ %in% c('brain','Brain_WM','Brain_GM')){
      testing_data <- cbind(dependent_variable = data_matrix_test[, cc], x_test_data[,c('eTIV',SiteCols)])
    }else{
      testing_data <- cbind(dependent_variable = data_matrix_test[, cc], x_test_data[,c('Height',SiteCols)])
    }
    # Predict the effect of covariates on the testing data using the training model
    predicted_values_test <- predict(best_covModel[[cc]], newdata = testing_data)
    # Calculate residuals for the testing data
    residuals_matrix_test[,cc] <- testing_data$dependent_variable - predicted_values_test
  }
  df[,unlist(sapply(colnames(df),function(x) grepl('^X',x)), use.names = F)] <- residuals_matrix_test
  df <- df[, !(names(df) %in% c(unlist(sapply(colnames(df),function(x) {if(grepl('Site',x)) x}), use.names = F),'eTIV','Height','Sex'))]
  
  
  # normalization and standardization
  df[,2:ncol(df)-1] <- predict(best_scaler, df[,2:ncol(df)-1])



  X_all <- df[, -c(1, ncol(df))]
  y_all_pred <- as.numeric(predict(best_model, newx = as.matrix(X_all), s = best_lambda))
  y_all_pred_corrected <- (y_all_pred - coef(best_correction_model)[1]) / coef(best_correction_model)[2]
  all_delta_corrected <- y_all_pred_corrected - df$Age
  df$pre_age <- y_all_pred
  df$pre_age_corrected <- y_all_pred_corrected
  df$delta_corrected <- all_delta_corrected
  df_all <- rbind(df_HC, df)
  write_csv(data.frame(df_all), paste0('/data/ageModeling/', organ, '_age_prediction_all.csv'))
}
sink()


# annotate the field name of weights
organs <- c('brain', 'Brain_WM', 'Brain_GM', 'heart', 'kidney', 'liver', 'Pancreas', 'OCT', 'sub','bone', 'boneSub','body','overall')
fieldnames <- read.table("/data/all_traits_ID_name.csv", sep=',', header=T)
fieldnames <- fieldnames[,c('Traits', 'Desc')]
colnames(fieldnames) <- c('field','Desc') 
fieldnames$field <- sapply(fieldnames$field, function(x) paste0('X',gsub('-','.',x)))
for(gender in genders){
  for(organ in organs){
    lasso_weights <- read.table(paste0('/data/ageModeling/', organ, 'Age_LASSO_weights','.csv'),sep=',',header=T)
    output <- merge(lasso_weights, fieldnames,by='field',all.x=T)
    write.table(output, paste0('/data/ageModeling/', organ, 'Age_LASSO_weights','_annotated.csv'), sep=',', row.names=F)
  }
}
