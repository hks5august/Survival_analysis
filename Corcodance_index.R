library(survival)
library(survminer)
library(dplyr)
library(survcomp)

set.seed(7)

train <- read.csv(file='final_train_pred_results',header =TRUE, sep = ",", dec = ".", row.names = 1);
################## remove samples where OS.time or DFI time is NA ##########
train<-subset(train,DFS.time!="nan")


test <- read.csv(file='final_test_pred_results',header =TRUE, sep = ",", dec = ".", row.names = 1);
################## remove samples where OS.time or DFI time is NA ##########
test<-subset(test,DFS.time!="nan")



surv_train <- Surv(train$DFS.time, train$DFS)
surv_test <- Surv(test$DFS.time, test$DFS)
#fit <- coxph(surv ~ group, data= sample.data)
#coxPredict <- predict(fit, data=sample.data, type="risk")  
predict_tr<-train$pred
predict_te<-test$pred

tr_CI<-concordance.index(x=predict_tr, surv.time=train$DFS.time, surv.event=train$DFS, method="noether")
#### x=predictions
tr_CI$c.index

te_CI<-concordance.index(x=predict_te, surv.time=test$DFS.time, surv.event=test$DFS, method="noether")

te_CI$c.index
