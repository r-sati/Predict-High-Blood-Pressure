library(data.table)
library(plyr)
library(dplyr)
library(caTools)
library(mlbench)
library(caret)
library(Boruta)
library(ROCR)

# Read input data 
inp_diet = read.csv("/Users/rsati/Documents/Personal/GaTech/MGT 6203/Project/Team-68-main/Data/Cleaned/diet_clean.csv")
inp_exam = read.csv("/Users/rsati/Documents/Personal/GaTech/MGT 6203/Project/Team-68-main/Data/Cleaned/examination_cleaned.csv")
inp_labs= read.csv("/Users/rsati/Documents/Personal/GaTech/MGT 6203/Project/Team-68-main/Data/Cleaned/labs_clean.csv")
inp_med= read.csv("/Users/rsati/Documents/Personal/GaTech/MGT 6203/Project/Team-68-main/Data/Cleaned/medications_cleaned.csv")
inp_target= read.csv("/Users/rsati/Documents/Personal/GaTech/MGT 6203/Project/Team-68-main/Data/Cleaned/mgt_6203_labeled_bp_data.csv")

inp_demo_seq= read.csv("/Users/rsati/Documents/Personal/GaTech/MGT 6203/Project/Team-68-main/Data/Cleaned/non_null_demographic_seqn.csv")
inp_demo_orig = read.csv("/Users/rsati/Documents/Personal/GaTech/MGT 6203/Project/Team-68-main/Data/Original/demographic.csv")

inp_demo<- inner_join(inp_demo_orig, inp_demo_seq, by = "SEQN") 


inp_all<- inner_join(inp_diet, inp_exam, by = "SEQN") 
inp_all<- inner_join(inp_all, inp_labs, by = "SEQN") 
inp_all<- inner_join(inp_all, inp_demo, by = "SEQN") 

# clean data 
# get percent of null values & filter for >= 75% data points
nm<-sapply(inp_all, function(x) sum(is.na(x)))/nrow(inp_all)
inp_all<- inp_all[,which(nm<0.25)]

# getting median of each column 
col_median <- apply(inp_all, 2, median, na.rm=TRUE)

# impute missing Values with median of Column
for(i in colnames(inp_all)){
  inp_all[,i][is.na(inp_all[,i])] <- col_median[i]
}

# Train & Test data split. 70% train and 30% test.

sample <- sample.split(inp_all$elevated_bp_flag, SplitRatio = 0.7)
train <- subset(inp_all, sample == TRUE)
target <- train$elevated_bp_flag
train <- subset(train, select = -c(elevated_bp_flag) )
test  <- subset(inp_all, sample == FALSE)
test_target <- test$elevated_bp_flag
test <- subset(test, select = -c(elevated_bp_flag) )


# Features Selection using boruta algorithm. 

set.seed(111)
#boruta <- Boruta(target ~ ., data = train, doTrace = 2, maxRuns = 500)
print(boruta)

plot(boruta, las = 2, cex.axis = 0.7)

final.boruta <- TentativeRoughFix(boruta)
getSelectedAttributes(final.boruta, withTentative = F)
boruta.df <- attStats(final.boruta)

# filter for important features. Out of 151 variables, 70 are selected as important.
selected_columns <- rownames(boruta.df[boruta.df$decision == "Confirmed",])
train_new = train[selected_columns]
test_new = test[selected_columns]


# Logistic regression model 

lm <- glm(target ~.,family=binomial(link='logit'),data=train_new)
summary(lm)


# model fitting
fitted.results <- predict(lm,newdata=test,type='response')
fitted.results <- ifelse(fitted.results > 0.5,1,0)

# accuracy 
misClasificError <- mean(fitted.results != test_target)
print(paste('Accuracy',1-misClasificError))
# accuracy = 0.76762


# plot the ROC curve and calculate the AUC (area under the curve) 
p <- predict(lm, newdata=test, type="response")
pr <- prediction(p, test_target)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

#auc = 0.8414795

# write.csv(test,"/Users/rsati/Documents/Personal/GaTech/MGT 6203/Project/test_output.csv",row.names = FALSE)


### Random forest

install.packages("randomForest")
library(randomForest)

# build model
model_base <- randomForest(target ~ ., data = train_new, importance = TRUE)
model_base

# prediction and confusion matrix
pred_base <- predict(model_base, test_new,  type = "class")
pred_base <- ifelse(pred_base > 0.5,1,0)

confusionMatrix(as.factor(pred_base), as.factor(test_target))

# accuracy
mean(pred_base == test_target)                    

# plot the ROC curve and calculate the AUC (area under the curve) 
p <- predict(model_base, test_new)
pr <- prediction(p, test_target)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc


# Tune parameters

tune<- tuneRF(
  x = train_new,
  y = target,
  ntreeTry = 500,    # tuning number of trees
  stepFactor = 0.5,
  improve = 0.01,    # it is the relative improvement in OOB error and so the improvement of the RM depends on the 0.05
  trace = TRUE
)

tune

model_tuned <- randomForest(target ~ ., data = train_new, ntree = 500, mtry = 23,  importance = TRUE)
model_tuned

pred_tuned <- predict(model_tuned, test_new,  type = "class")
pred_tuned <- ifelse(pred_tuned > 0.5,1,0)

confusionMatrix(as.factor(pred_tuned), as.factor(test_target))

mean(pred_tuned == test_target)

table(pred_tuned,test_target)

# plot the ROC curve and calculate the AUC (area under the curve) 
p <- predict(model_tuned, test_new)
pr <- prediction(p, test_target)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

## important variable
varImpPlot(model_tuned)
importance(model_tuned)

## XG Boost 

install.packages("xgboost")
library(xgboost)

# Base XGBoost model
set.seed(20)
params <- list(booster = "gbtree", 
               objective = "binary:logistic")

xgb_base <- xgb.train (params = params,
                       data = train_dmatrix,
                       nrounds =1000,
                       print_every_n = 10,
                       eval_metric = "auc",
                       eval_metric = "error",
                       early_stopping_rounds = 50,
                       watchlist = list(train= train_dmatrix, val= valid_dmatrix))

# predictions and confusion matrix
predict_xgb_base <- predict(xgb_base, as.matrix(test_new))
fitted.results_xgb_base <- ifelse(predict_xgb_base > 0.5,1,0)

misClasificError_xgb_base <- mean(fitted.results_xgb_base != test_target)
print(paste('Accuracy',1-misClasificError_xgb_base))
# accuracy = 0.7676

confusionMatrix(as.factor(fitted.results_xgb_base), as.factor(test_target))

# plot the ROC curve and calculate the AUC (area under the curve) 
p <- predict(xgb_base, as.matrix(test_new))
pr <- prediction(p, test_target)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
#auc = 0.836051


# Create 1,000 rows with random hyperparameters

# Create empty lists
lowest_error_list = list()
parameters_list = list()

set.seed(20)
for (iter in 1:1000){
  param <- list(booster = "gbtree",
                objective = "binary:logistic",
                max_depth = sample(3:10, 1),
                eta = runif(1, .01, .3),
                subsample = runif(1, .7, 1),
                colsample_bytree = runif(1, .6, 1),
                min_child_weight = sample(0:10, 1)
  )
  parameters <- as.data.frame(param)
  parameters_list[[iter]] <- parameters
}

# Create object that contains all randomly created hyperparameters
parameters_df = do.call(rbind, parameters_list)


# Use randomly created parameters to create 1,000 XGBoost-models
for (row in 1:nrow(parameters_df)){
  set.seed(20)
  mdcv <- xgb.train(data=train_dmatrix,
                    booster = "gbtree",
                    objective = "binary:logistic",
                    max_depth = parameters_df$max_depth[row],
                    eta = parameters_df$eta[row],
                    subsample = parameters_df$subsample[row],
                    colsample_bytree = parameters_df$colsample_bytree[row],
                    min_child_weight = parameters_df$min_child_weight[row],
                    nrounds= 300,
                    eval_metric = "error",
                    early_stopping_rounds= 30,
                    print_every_n = 100,
                    watchlist = list(train= train_dmatrix, val= valid_dmatrix)
  )
  lowest_error <- as.data.frame(1 - min(mdcv$evaluation_log$val_error))
  lowest_error_list[[row]] <- lowest_error
}


# Create object that contains all accuracy's
lowest_error_df = do.call(rbind, lowest_error_list)

# Bind columns of accuracy values and random hyperparameter values
randomsearch = cbind(lowest_error_df, parameters_df)

# highest accuracy
max(randomsearch$`1 - min(mdcv$evaluation_log$val_error)`)

randomsearch_tuned

randomsearch_tuned[1,]$max_depth
randomsearch_tuned[1,]$eta
randomsearch_tuned[1,]$subsample
randomsearch_tuned[1,]$colsample_bytree
randomsearch_tuned[1,]$min_child_weight

# Tuned XGBoost model

set.seed(20)
params <- list(booster = "gbtree", 
               objective = "binary:logistic",
               max_depth = randomsearch_tuned[1,]$max_depth,
               eta = randomsearch_tuned[1,]$eta,
               subsample = randomsearch_tuned[1,]$subsample,
               colsample_bytree = randomsearch_tuned[1,]$colsample_bytree,
               min_child_weight = randomsearch_tuned[1,]$min_child_weight)

xgb_tuned <- xgb.train(params = params,
                       data = train_dmatrix,
                       nrounds =1000,
                       print_every_n = 10,
                       eval_metric = "auc",
                       eval_metric = "error",
                       early_stopping_rounds = 30,
                       watchlist = list(train= train_dmatrix, val= valid_dmatrix))


# Make prediction

# predictions and confusion matrix
predict_xgb_tuned <- predict(xgb_tuned, as.matrix(test_new))
fitted.results_xgb_tuned <- ifelse(predict_xgb_tuned > 0.5,1,0)

misClasificError_xgb_tuned <- mean(fitted.results_xgb_tuned != test_target)
print(paste('Accuracy',1-misClasificError_xgb_tuned))
# accuracy = 0.7676

confusionMatrix(as.factor(fitted.results_xgb_tuned), as.factor(test_target))

# plot the ROC curve and calculate the AUC (area under the curve) 
p <- predict(xgb_tuned, as.matrix(test_new))
pr <- prediction(p, test_target)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc
#auc = 0.8430882


# feature importances
importance_matrix = xgb.importance(colnames(train_dmatrix), model = xgb_tuned)
head(importance_matrix)

xgb.plot.importance(importance_matrix[1:10,])


