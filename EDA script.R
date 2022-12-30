library(data.table)

# Read input data 
inp_dt = read.csv("/Users/rsati/Documents/Personal/GaTech/MGT 6203/Project/labs.csv")

# get percent of null values & filter for >= 75% data points
nm<-sapply(inp_dt, function(x) sum(is.na(x)))/nrow(inp_dt)
dt<- inp_dt[,which(nm<0.25)]

summary(dt)
str(dt)


# getting median of each column 
col_median <- apply(dt, 2, median, na.rm=TRUE)

# impute missing Values with median of Column
for(i in colnames(dt)){
  dt[,i][is.na(dt[,i])] <- col_median[i]
}


# Remove Correlated Variables

cor_matrix <- cor(dt)                      # Correlation matrix
cor_matrix

cor_matrix_rm <- cor_matrix                  # Modify correlation matrix
cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0
diag(cor_matrix_rm) <- 0
cor_matrix_rm

dt_new <- dt[ , !apply(cor_matrix_rm,    # Remove highly correlated variables
                           2,
                           function(x) any(x > 0.80))]
ncol(dt_new)  

write.csv(dt_new,"/Users/rsati/Documents/Personal/GaTech/MGT 6203/Project/labs_clean.csv",row.names = FALSE)

