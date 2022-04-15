library(survival)
library(survminer)
library(dplyr)
library(survcomp)

library(survAUC)


########## cross validation ###########
set.seed(7)
train <- read.csv(file='final_train_predictions',header =TRUE, sep = ",", dec = ".", );
################## remove samples where OS.time or DFI time is NA ##########
train<-subset(train,DFS.time!="nan")

test <- read.csv(file="final_test_predictions",header =TRUE, sep = ",", dec = ".", );
################## remove samples where OS.time or DFI time is NA ##########
test<-subset(test,DFS.time!="nan")



ext_test <- read.csv(file="final_ext", header =TRUE, sep = ",", dec = ".", );
################## remove samples where OS.time or DFI time is NA ##########

ext_test<-subset(ext_test,DFS.time!="nan")


tr_CI<-concordance.index(x=train$Ridge.Classifier, surv.time=train$DFS.time, surv.event=train$DFS, method="noether")
#### x=predictions

tr_CI

D_index_tr<-D.index(x=train$pred1, surv.time=train$DFS.time, surv.event=train$DFS, alpha = 0.05,
                    method.test = c("logrank"), na.rm = T)
D_index_tr

HR_tr<-hazard.ratio(x=train$Ridge.Classifier, surv.time=train$DFS.time, surv.event=train$DFS, alpha = 0.05,
                    method.test = c("logrank"), na.rm = T)

HR_tr



########### Test ##############

te_CI<-concordance.index(x=test$Ridge.Classifier, surv.time=test$DFS.time, surv.event=test$DFS, method="noether")
#### x=predictions

te_CI

D_index_te<-D.index(x=test$svc.rbf_w, surv.time=test$DFS.time, surv.event=test$DFS, alpha = 0.05,
                    method.test = c("logrank"), na.rm = T)
D_index_te

HR_te<-hazard.ratio(x=test$Ridge.Classifier, surv.time=test$DFS.time, surv.event=test$DFS, alpha = 0.05,
                    method.test = c("logrank"), na.rm = T)

HR_te



########### Test ##############

ext_te_CI<-concordance.index(x=ext_test$pred, surv.time=ext_test$DFS.time, surv.event=ext_test$DFS, method="noether")
#### x=predictions

ext_te_CI
#ext_te_CI$c.index
#ext_te_CI$p.value

ext_D_index_te<-D.index(x=ext_test$svc.rbf_w, surv.time=ext_test$DFS.time, surv.event=ext_test$DFS, alpha = 0.05,
                    method.test = c("logrank"), na.rm = T)
ext_D_index_te

ext_HR_te<-hazard.ratio(x=ext_test$pred, surv.time=ext_test$DFS.time, surv.event=ext_test$DFS, alpha = 0.05,
                    method.test = c("logrank"), na.rm = T)

ext_HR_te



######## Brier Score ####

dd_tr <- data.frame("time"=train$DFS.time, "event"=train$DFS, "score"=train$Ridge.Classifier)
dd_te <- data.frame("time"=test$DFS.time, "event"=test$DFS, "score"=test$Ridge.Classifier)
dd_ext <- data.frame("time"=ext_test$DFS.time, "event"=ext_test$DFS, "score"=ext_test$pred)


tr_Brier_score <- sbrier.score2proba(data.tr=dd_tr, data.ts=dd_tr, method="cox")
te_Brier_score <- sbrier.score2proba(data.tr=dd_te, data.ts=dd_te, method="cox")
ext_te_Brier_score <- sbrier.score2proba(data.tr=dd_ext, data.ts=dd_ext, method="cox")

tr_Brier_score$bsc.integrated
te_Brier_score$bsc.integrated
ext_te_Brier_score$bsc.integrated

