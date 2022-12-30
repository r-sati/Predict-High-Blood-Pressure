library(data.table)
library(plyr)
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




