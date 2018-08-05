####################################
# MISSING DATA
####################################

#clean environment 
rm(list=ls())


#set working directory
setwd("~/Desktop/summer_projects/model_based_missing_data/data")


#load libraries
library(MASS)
library(MVN)
library(matrixcalc)
library(tidyr)
library(mice)
library(Amelia)
library(Zelig)
library(broom)
library(tibble)


###########################################################
#simulate complte synthetic multivariate normal dataset
###########################################################
#correlation matrix
set.seed(1)
sig <- matrix(runif(121, 0, 1000), ncol=11) 
sig <- t(sig) %*% sig

colnames(sig) = c("Y", "X1", "X2", "X3" ,"X4", "X5", "X6", "X7", "X8", "X9", "X10")
rownames(sig) = colnames(sig)


# vector of means
set.seed(1)
mu<-c(round(runif(11, 0,1000),3))

set.seed(1)
sim_comp<-as.data.frame(mvrnorm(n = 500, mu = mu, Sigma = sig))

write.csv(sim_comp, "sim_comp.csv", row.names = F)
sim_comp<-read.csv("sim_comp.csv")

###############################################################################
# Create datasets with specified percent of observations with missing data    #
###############################################################################
# Randomly delete values from 10 percent of rows 
# Section has been commented out to ensure same dataset is used everytime

# set.seed(1) 
# rows<-sample(1:nrow(sim_comp), 0.10*nrow(sim_comp), replace = FALSE) #sample indices 10 percent of rows

# sim_miss_10<-sim_comp #make copy of dataset 

#for each row sample 1 to 6 columns to make NA
# for(row in rows){
#   cols<-sample(1:ncol(sim_miss_10), sample(1:6,1), replace = F) 
#   sim_miss_10[row, cols]<-NA
# }

#write sim_miss_10 as CVS
# write.csv(sim_miss_10, "sim_miss_10.csv", row.names = F)

#load sim_miss_10
sim_miss_10<-read.csv("sim_miss_10.csv")


# Randomly delete values from 25 percent of rows 
# Section has been commented out to ensure same dataset is used everytime

# set.seed(1)
#  rows<-sample(1:nrow(sim_comp), 0.25*nrow(sim_comp), replace = FALSE) #sample indices 25 percent of rows

# sim_miss_25<-sim_comp #make copy of dataset

#for each row sample 1 to 6 columns to make NA
# for(row in rows){
#   cols<-sample(1:ncol(sim_miss_25), sample(1:6,1), replace = F)
#   sim_miss_25[row, cols]<-NA
# }

#write out dataset
# write.csv(sim_miss_25, "sim_miss_25.csv", row.names = F)

#load sim_miss_25
sim_miss_25<-read.csv("sim_miss_25.csv")


# Randomly delete values from 50 percent of rows 
# Section has been commented out to ensure same dataset is used everytime

# set.seed(3)
# rows<-sample(1:nrow(sim_comp), 0.5*nrow(sim_comp), replace = FALSE)

# sim_miss_50<-sim_comp

# for(row in rows){
#   cols<-sample(1:ncol(sim_miss_50), sample(1:6,1), replace = F)  #sample indices 50 percent of rows
#   sim_miss_50[row, cols]<-NA
# }


#write out dataset
# write.csv(sim_miss_50, "sim_miss_50.csv", row.names = F)

#read in dataset
sim_miss_50<-read.csv("sim_miss_50.csv")


################################################## 
# Complete case Analysis: Means
################################################## 
#get mean of each column for each dataset
colMeans(sim_comp)
colMeans(sim_miss_10, na.rm = T)
colMeans(sim_miss_25, na.rm = T)
colMeans(sim_miss_50, na.rm = T)


################################################## 
# Complete case Analysis: standard deviations
################################################## 
#get standard deviation of each column for each dataset
apply(sim_comp, 2, sd)
apply(sim_miss_10[complete.cases(sim_miss_10),], 2, sd)
apply(sim_miss_25[complete.cases(sim_miss_25),], 2, sd)
apply(sim_miss_50[complete.cases(sim_miss_50),], 2, sd)


################################################## 
# Complete case Analysis: Regression 3 predictors
################################################## 
#run regression with X1, X2 and X4 as predictiors for each dataset
lm_comp<-lm(Y~ X1 + X2 + X4, data = sim_comp)
lm_miss_10<-lm(Y~ X1 + X2 + X4, data = sim_miss_10)
lm_miss_25<-lm(Y~ X1 + X2 + X4, data = sim_miss_25)
lm_miss_50<-lm(Y~ X1 + X2 + X4, data = sim_miss_50)


#reformat results
lm_miss_10 <- tidy(lm_miss_10)
lm_miss_25 <- tidy(lm_miss_25)
lm_miss_50 <- tidy(lm_miss_50)


#write out results in csv
write.csv(lm_miss_10,"lm_miss_10.csv", row.names = F)
write.csv(lm_miss_25,"lm_miss_25.csv", row.names = F)
write.csv(lm_miss_50,"lm_miss_50.csv", row.names = F)



####################################################
# Multiple imputation: EM
####################################################
#run algorithm missing 10 percent 
ptm <- proc.time()  #start time 
amelia_10 <- amelia(sim_miss_10, m=5) #algorithm 
duration_amelia_10<-proc.time() - ptm  #end time

ptm <- proc.time() #start time 
amelia_25 <- amelia(sim_miss_25, m=5) #algorithm 
duration_amelia_25<-proc.time() - ptm #end time

ptm <- proc.time()  #start time 
amelia_50 <- amelia(sim_miss_50, m=5) #algorithm 
duration_amelia_50<-proc.time() - ptm #end time


duration_amelia_10<-t(data.matrix(duration_amelia_10))
duration_amelia_25<-t(data.matrix(duration_amelia_25))
duration_amelia_50<-t(data.matrix(duration_amelia_50))

#write out duration for each level of missingness 
write.csv(duration_amelia_10, "duration_amelia_10.csv", row.names = F)
write.csv(duration_amelia_25, "duration_amelia_25.csv", row.names = F)
write.csv(duration_amelia_50, "duration_amelia_50.csv", row.names = F)



#calculate column means for each level of missingness over all imputed datasets 
amelia_10_colMeans<-rowMeans(sapply(amelia_10$imputations, colMeans))
amelia_25_colMeans<-rowMeans(sapply(amelia_25$imputations, colMeans))
amelia_50_colMeans<-rowMeans(sapply(amelia_50$imputations, colMeans))

#write out column means
write.csv(amelia_10_colMeans, "amelia_10_colMeans.csv", row.names = F)
write.csv(amelia_25_colMeans, "amelia_25_colMeans.csv", row.names = F)
write.csv(amelia_50_colMeans, "amelia_50_colMeans.csv", row.names = F)

#calculate column standard deviations for each level of missingness over all imputed datasets 
amelia_10_sd<-rowMeans(sapply(amelia_10$imputations, sapply, sd))
amelia_25_sd<-rowMeans(sapply(amelia_25$imputations, sapply, sd))
amelia_50_sd<-rowMeans(sapply(amelia_50$imputations, sapply, sd))

#write out column standard deviations
write.csv(amelia_10_sd, "amelia_10_sd.csv", row.names = F)
write.csv(amelia_10_sd, "amelia_10_sd.csv", row.names = F)
write.csv(amelia_10_sd, "amelia_10_sd.csv", row.names = F)


#run regressions on imputed datasets for each leve of missingness
lm_amelia_10 <- zelig(Y ~ X1 + X2 + X4, data = amelia_10, model = "ls", cite=F)
lm_amelia_25 <- zelig(Y ~ X1 + X2 + X4, data = amelia_25, model = "ls", cite=F)
lm_amelia_50 <- zelig(Y ~ X1 + X2 + X4, data = amelia_50, model = "ls", cite=F)

#function to extract coefficients standard errors and pvalues from amelia object
extract_info<-function(amelia_object){
  amelia_coef<-rowMeans(sapply(amelia_object$get_coef(),sapply, mean))
  amelia_se<-rowMeans(sapply(amelia_object$get_se(),sapply, mean))
  amelia_pval<-rowMeans(sapply(amelia_object$get_pvalue(),sapply, mean))
  
  combined<-t(rbind(amelia_coef,amelia_se,amelia_pval))
  colnames(combined)<-c("Estimate", "Std.Error", "Pr(>|z|)")
  
  return(combined)
}

#extract information from the three levels of missingness
lm_amelia_10<-extract_info(lm_amelia_10)
lm_amelia_25<-extract_info(lm_amelia_25)
lm_amelia_50<-extract_info(lm_amelia_50)

#write out regression results
write.csv(lm_amelia_10, "lm_amelia_10.csv")
write.csv(lm_amelia_25, "lm_amelia_25.csv")
write.csv(lm_amelia_50, "lm_amelia_50.csv")



####################################################
# Multiple imputation: PMM
####################################################
#select all columns except X5
sim_miss_10<-sim_miss_10[,c("Y","X1","X2","X4", "X6",
                            "X7", "X8", "X9","X10")] #if all included warnings produced


#select all columns except X5
sim_miss_25<-sim_miss_25[,c("Y","X1","X2","X4", "X6",
                            "X7", "X8", "X9", "X10")] #if all included warnings produced

#select all columns except X5
sim_miss_50<-sim_miss_50[,c("Y","X1","X2","X4", "X6",
                            "X7", "X8", "X9", "X10")] #if all included warnings produced


#impute values 5 times with a maximum of 50 iterations using predicted mean matching
ptm <- proc.time() #start time
pmm_10<-mice(sim_miss_10, m=5, maxit = 50, method = 'pmm', seed = 500)
duration_pmm_10<-proc.time() - ptm # Stop time 

ptm <- proc.time() #start time
pmm_25<-mice(sim_miss_25, m=5, maxit = 50, method = 'pmm', seed = 500)
duration_pmm_25<-proc.time() - ptm # Stop time 

ptm <- proc.time() #start time
pmm_50<-mice(sim_miss_50, m=5, maxit = 50, method = 'pmm', seed = 500)
duration_pmm_50<-proc.time() - ptm # Stop time 


duration_pmm_10<-t(data.matrix(duration_pmm_10))
duration_pmm_25<-t(data.matrix(duration_pmm_25))
duration_pmm_50<-t(data.matrix(duration_pmm_50))


#write out duration
write.csv(duration_pmm_10, "duration_pmm_10.csv", row.names = F)
write.csv(duration_pmm_25, "duration_pmm_25.csv", row.names = F)
write.csv(duration_pmm_50, "duration_pmm_50.csv", row.names = F)


#extract column means from imputed data 
pmm_10_colMeans<-as.data.frame(colMeans(mice::complete(pmm_10, 'long')[,c(3:11)]))
pmm_10_colMeans<-as.data.frame(cbind(Variable = row.names(pmm_10_colMeans), pmm_10_colMeans))
colnames(pmm_10_colMeans)<-c("Variable", "Column Mean")

pmm_25_colMeans<-as.data.frame(colMeans(mice::complete(pmm_25, 'long')[,c(3:11)]))
pmm_25_colMeans<-as.data.frame(cbind(Variable = row.names(pmm_25_colMeans), pmm_25_colMeans))
colnames(pmm_25_colMeans)<-c("Variable", "Column Mean")

pmm_50_colMeans<-as.data.frame(colMeans(mice::complete(pmm_50, 'long')[,c(3:11)]))
pmm_50_colMeans<-as.data.frame(cbind(Variable = row.names(pmm_50_colMeans), pmm_50_colMeans))
colnames(pmm_50_colMeans)<-c("Variable", "Column Mean")


#write out aggregate ColMeans for all 5 complete datasets
write.csv(pmm_10_colMeans, "pmm_10_colMeans.csv", row.names = F)
write.csv(pmm_25_colMeans, "pmm_25_colMeans.csv", row.names = F)
write.csv(pmm_50_colMeans, "pmm_50_colMeans.csv", row.names = F)


#extract standard deviations from imputed data 
pmm_10_colSD<-sapply(mice::complete(pmm_10, 'long')[,c(3:11)], sd)
pmm_25_colSD<-sapply(mice::complete(pmm_25, 'long')[,c(3:11)], sd)
pmm_50_colSD<-sapply(mice::complete(pmm_50, 'long')[,c(3:11)], sd)

#write out column SD
write.csv(pmm_10_colSD,"pmm_10_colSD.csv", row.names = F)
write.csv(pmm_25_colSD,"pmm_25_colSD.csv", row.names = F)
write.csv(pmm_50_colSD,"pmm_50_colSD.csv", row.names = F)



#with imputed datasets run regression with X1, X2, X4 as predictors
lm_pmm_10<-with(pmm_10, lm(Y~X1 + X2 + X4))

#with imputed datasets run regression with X1, X2, X4 as predictors
lm_pmm_25<-with(pmm_25, lm(Y~X1 + X2 + X4))

#with imputed datasets run regression with X1, X2, X4 as predictors
lm_pmm_50<-with(pmm_50, lm(Y~X1 + X2 + X4))


#save results of regressions on imputed data as dataframes
lm_pmm_10<-as.data.frame(summary(pool(lm_pmm_10)))
lm_pmm_25<-as.data.frame(summary(pool(lm_pmm_25)))
lm_pmm_50<-as.data.frame(summary(pool(lm_pmm_50)))


#write out results of regressions on inputed data as CSV
write.csv(lm_pmm_10,"lm_pmm_10.csv")
write.csv(lm_pmm_25,"lm_pmm_25.csv")
write.csv(lm_pmm_50,"lm_pmm_50.csv")




###########################################
# Read data set
###########################################
rm(list=ls())

#load in cars data
cars_comp<-read.table("auto-mpg.data",header = FALSE)

#add column names
names(cars_comp)<-c("mpg","cylinders","displacement","horsepower","weight",
               "acceleration", "model_year", "origin", "car_name")

#select numeric columns
cars_comp<-cars_comp[,c(1:7)]

#convert columns from character to numeric
cars_comp<-sapply(cars_comp, as.numeric)

#convert from matrix to dataframe
cars_comp<-as.data.frame(cars_comp)

#create histograms to show data not exactly multivariate normal
hist(cars_comp$mpg)
hist(cars_comp$cylinders)
hist(cars_comp$displacement)
hist(cars_comp$horsepower)
hist(cars_comp$weight)
hist(cars_comp$acceleration)
hist(cars_comp$model_year)


###############################################################################
# Create datasets with specified percent of observations with missing data    #
###############################################################################
# Randomly delete values from 10 percent of rows 
# Section has been commented out to ensure same dataset is used everytime

# set.seed(1)
# rows<-sample(1:nrow(cars_comp), 0.10*nrow(cars_comp), replace = FALSE)

# cars_miss_10<-cars_comp

# for(row in rows){
#   cols<-sample(1:ncol(cars_miss_10), sample(1:6,1), replace = F)
#   cars_miss_10[row, cols]<-NA
# }

# write.csv(cars_miss_10, "cars_miss_10.csv", row.names = F)

cars_miss_10<-read.csv("cars_miss_10.csv")



# Randomly delete values from 25 percent of rows 
# Section has been commented out to ensure same dataset is used everytime

# set.seed(2)
# rows<-sample(1:nrow(cars_comp), 0.25*nrow(cars_comp), replace = FALSE)

# cars_miss_25<-cars_comp

# for(row in rows){
#   cols<-sample(1:ncol(cars_miss_25), sample(1:6,1), replace = F)
#   cars_miss_25[row, cols]<-NA
# }

# write.csv(cars_miss_25, "cars_miss_25.csv", row.names = F)

cars_miss_25<-read.csv("cars_miss_25.csv")


# Randomly delete values from 50 percent of rows 
# Section has been commented out to ensure same dataset is used everytime

# set.seed(3)
# rows<-sample(1:nrow(cars_comp), 0.5*nrow(cars_comp), replace = FALSE)

# cars_miss_50<-cars_comp

# for(row in rows){
#   cols<-sample(1:ncol(cars_miss_50), sample(1:6,1), replace = F)
#   cars_miss_50[row, cols]<-NA
# }

# write.csv(cars_miss_50, "cars_miss_50.csv", row.names = F)

cars_miss_50<-read.csv("cars_miss_50.csv")



################################################## 
# Complete case Analysis: Means
##################################################
#get mean of each column for each dataset
colMeans(cars_comp)
colMeans(cars_miss_10, na.rm = T)
colMeans(cars_miss_25, na.rm = T)
colMeans(cars_miss_50, na.rm = T)


################################################## 
# Complete case Analysis: standard deviations
################################################## 
#get standard deviation of each column for each dataset
apply(cars_comp, 2, sd)
apply(cars_miss_10[complete.cases(cars_miss_10),], 2, sd)
apply(cars_miss_25[complete.cases(cars_miss_25),], 2, sd)
apply(cars_miss_50[complete.cases(cars_miss_50),], 2, sd)


################################################## 
# Complete case Analysis: Regression all predictors
################################################## 
#run regression using all predictors for each dataset
lm_cars_comp<-lm(mpg~ . , data = cars_comp)
lm_cars_miss_10<-lm(mpg~ . , data = cars_miss_10)
lm_cars_miss_25<-lm(mpg~ . , data = cars_miss_25)
lm_cars_miss_50<-lm(mpg~ . , data = cars_miss_50)




####################################################
# Multiple imputation: EM
####################################################

#run algorithm
cars_amelia_10 <- amelia(cars_miss_10, m=5) #algorithm 
cars_amelia_25 <- amelia(cars_miss_25, m=5) #algorithm 
cars_amelia_50 <- amelia(cars_miss_50, m=5) #algorithm 



#calculate column means for each level of missingness over all imputed datasets 
cars_amelia_10_colMeans<-rowMeans(sapply(cars_amelia_10$imputations, colMeans))
cars_amelia_25_colMeans<-rowMeans(sapply(cars_amelia_25$imputations, colMeans))
cars_amelia_50_colMeans<-rowMeans(sapply(cars_amelia_50$imputations, colMeans))


#write out column means
write.csv(cars_amelia_10_colMeans, "cars_amelia_10_colMeans.csv", row.names = F)
write.csv(cars_amelia_25_colMeans, "cars_amelia_25_colMeans.csv", row.names = F)
write.csv(cars_amelia_50_colMeans, "cars_amelia_50_colMeans.csv", row.names = F)


#calculate column standard deviations for each level of missingness over all imputed datasets 
cars_amelia_10_sd<-rowMeans(sapply(cars_amelia_10$imputations, sapply, sd))
cars_amelia_25_sd<-rowMeans(sapply(cars_amelia_25$imputations, sapply, sd))
cars_amelia_50_sd<-rowMeans(sapply(cars_amelia_50$imputations, sapply, sd))

#write out column standard deviations
write.csv(cars_amelia_10_sd, "cars_amelia_10_sd.csv", row.names = F)
write.csv(cars_amelia_10_sd, "cars_amelia_10_sd.csv", row.names = F)
write.csv(cars_amelia_10_sd, "cars_amelia_10_sd.csv", row.names = F)


#run regressions on imputed datasets for each leve of missingness
cars_lm_amelia_10 <- zelig(mpg ~ cylinders + displacement + horsepower + weight+ acceleration+model_year, data = cars_amelia_10, model = "ls", cite=F)
cars_lm_amelia_25 <- zelig(mpg ~ cylinders + displacement + horsepower + weight+ acceleration+model_year, data = cars_amelia_25, model = "ls", cite=F)
cars_lm_amelia_50 <- zelig(mpg ~ cylinders + displacement + horsepower + weight+ acceleration+model_year, data = cars_amelia_50, model = "ls", cite=F)

#function to extract coefficients standard errors and pvalues from amelia object
extract_info<-function(amelia_object){
  amelia_coef<-rowMeans(sapply(amelia_object$get_coef(),sapply, mean))
  amelia_se<-rowMeans(sapply(amelia_object$get_se(),sapply, mean))
  amelia_pval<-rowMeans(sapply(amelia_object$get_pvalue(),sapply, mean))
  
  combined<-t(rbind(amelia_coef,amelia_se,amelia_pval))
  colnames(combined)<-c("Estimate", "Std.Error", "Pr(>|z|)")
  
  return(combined)
}

#extract information from the three levels of missingness
cars_lm_amelia_10<-extract_info(cars_lm_amelia_10)
cars_lm_amelia_25<-extract_info(cars_lm_amelia_25)
cars_lm_amelia_50<-extract_info(cars_lm_amelia_50)

#write out regression results
write.csv(cars_lm_amelia_10, "cars_lm_amelia_10.csv")
write.csv(cars_lm_amelia_25, "cars_lm_amelia_25.csv")
write.csv(cars_lm_amelia_50, "cars_lm_amelia_50.csv")




####################################################
# Multiple imputation: PMM
####################################################
#impute values 5 times with a maximum of 50 iterations using predicted mean matching
cars_pmm_10<-mice(cars_miss_10, m=5, maxit = 50, method = 'pmm', seed = 500)
cars_pmm_25<-mice(cars_miss_25, m=5, maxit = 50, method = 'pmm', seed = 500)
cars_pmm_50<-mice(cars_miss_50, m=5, maxit = 50, method = 'pmm', seed = 500)



#extract column means from each level of missingness
cars_pmm_10_colMeans<-as.data.frame(colMeans(mice::complete(cars_pmm_10, 'long')[,c(3:9)]))
names(cars_pmm_10_colMeans)<-c("Column Mean")

cars_pmm_25_colMeans<-as.data.frame(colMeans(mice::complete(cars_pmm_25, 'long')[,c(3:9)]))
names(cars_pmm_25_colMeans)<-c("Column Mean")

cars_pmm_50_colMeans<-as.data.frame(colMeans(mice::complete(cars_pmm_50, 'long')[,c(3:9)]))
names(cars_pmm_50_colMeans)<-c("Column Mean")

#write out column means
write.csv(cars_pmm_10_colMeans, "cars_pmm_10_colMeans.csv")
write.csv(cars_pmm_25_colMeans, "cars_pmm_25_colMeans.csv")
write.csv(cars_pmm_50_colMeans, "cars_pmm_50_colMeans.csv")

#extract standard deviations from imputed data 
cars_pmm_10_colSD<-sapply(mice::complete(cars_pmm_10, 'long')[,c(3:9)], sd)
cars_pmm_25_colSD<-sapply(mice::complete(cars_pmm_25, 'long')[,c(3:9)], sd)
cars_pmm_50_colSD<-sapply(mice::complete(cars_pmm_50, 'long')[,c(3:9)], sd)

#write out column SD
write.csv(cars_pmm_10_colSD,"pmm_10_colSD.csv", row.names = F)
write.csv(cars_pmm_25_colSD,"pmm_25_colSD.csv", row.names = F)
write.csv(cars_pmm_50_colSD,"pmm_50_colSD.csv", row.names = F)



#with imputed data run regression using all predictors
lm_cars_pmm_10<-with(cars_pmm_10, lm(mpg~cylinders + displacement + horsepower+
                                           weight+acceleration+model_year))

lm_cars_pmm_25<-with(cars_pmm_25, lm(mpg~cylinders + displacement + horsepower+
                                           weight+acceleration+model_year))

lm_cars_pmm_50<-with(cars_pmm_50, lm(mpg~cylinders + displacement + horsepower+
                                           weight+acceleration+model_year))



#save results of regressions on imputed data as dataframes
lm_cars_pmm_10<-as.data.frame(summary(pool(lm_cars_pmm_10)))
lm_cars_pmm_25<-as.data.frame(summary(pool(lm_cars_pmm_25)))
lm_cars_pmm_50<-as.data.frame(summary(pool(lm_cars_pmm_50)))


#write out results as CSV
write.csv(lm_cars_pmm_10,"lm_cars_pmm_10.csv")
write.csv(lm_cars_pmm_25,"lm_cars_pmm_25.csv")
write.csv(lm_cars_pmm_50,"lm_cars_pmm_50.csv")




