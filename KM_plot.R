library(survival)
library(survminer)
library(dplyr)
library(survcomp)
library("survminer")
library(survAUC)


########## Scripts for KM plot ###########
set.seed(7)
train <- read.csv(file='final_tr',header =TRUE, sep = ",", dec = ".", );
################## remove samples where OS.time or DFI time is NA ##########
train<-subset(train,DFS.time!="nan")

test <- read.csv(file="final_te",header =TRUE, sep = ",", dec = ".", );
################## remove samples where OS.time or DFI time is NA ##########
test<-subset(test,DFS.time!="nan")


ext_test <- read.csv(file="final_ext_te.csv", header =TRUE, sep = ",", dec = ".", );
################## remove samples where OS.time or DFI time is NA ##########

ext_test<-subset(ext_test,DFS.time1!="nan")




surv_train <- Surv(train$DFS.time, train$DFS)
#fit_tr <- coxph(surv_train ~ group, data= train)
#fit_tr <- coxph(surv_train ~ train$pred, data= train)

fit_tr <- survfit(surv_train ~ train$Risk, data= train)
fit_tr_age <- survfit(surv_train ~ (train$age <median(train$age)), data= train)
fit_tr_gender <- survfit(surv_train ~ train$gender, data= train)
fit_tr_stage <- survfit(surv_train ~ train$stage, data= train)
#fit_tr <- survfit(surv_train ~ train$pred, data= train)
#fit_tr <- coxph(surv_train ~ train$ETREES, data= train)
#fit_tr <- coxph(surv_train ~ train$svc.rbf_w, data= train)
fit_tr

median(train$age)

ggsurvplot(fit_tr, pval = TRUE, censor = TRUE, conf.int = TRUE, palette = c("blue", "red"), break.time.by = 365, risk.table = TRUE, risk.table.height = 0.5 ) #Useful when you have multiple groups 

#ggsurvplot(fit_tr_age)
ggsurvplot(fit_tr_age, pval = TRUE, censor = TRUE,  break.time.by = 365, risk.table = TRUE, risk.table.height = 0.5 ) #Useful when you have multiple groups 

ggsurvplot(fit_tr_gender, pval = TRUE, censor = TRUE,  break.time.by = 365,  risk.table = TRUE, risk.table.height = 0.5 ) #Useful when you have multiple groups 

ggsurvplot(fit_tr_stage, pval = TRUE, censor = TRUE, conf.int = F,  break.time.by = 365, risk.table = TRUE, risk.table.height = 0.5 ) #Useful when you have multiple groups 

######## Test Data ###########
surv_test <- Surv(test$DFS.time, test$DFS)
fit_te <- survfit(surv_test ~ test$Risk, data= test)

fit_te_age <- survfit(surv_test ~ (test$age <median(test$age)), data= test)
fit_te_gender <- survfit(surv_test ~ test$gender, data= test)
fit_te_stage <- survfit(surv_test ~ test$stage, data= test)

ggsurvplot(fit_te, pval = TRUE, censor = TRUE, conf.int = TRUE, palette = c("blue", "red"), break.time.by = 365, risk.table = TRUE, risk.table.height = 0.5 ) #Useful when you have multiple groups 

#ggsurvplot(fit_tr_age)
ggsurvplot(fit_te_age, pval = TRUE, censor = TRUE,  break.time.by = 365, risk.table = TRUE, risk.table.height = 0.5 ) #Useful when you have multiple groups 

ggsurvplot(fit_te_gender, pval = TRUE, censor = TRUE,  break.time.by = 365,  risk.table = TRUE, risk.table.height = 0.5 ) #Useful when you have multiple groups 

ggsurvplot(fit_te_stage, pval = TRUE, censor = TRUE, conf.int = F,  break.time.by = 365, risk.table = TRUE, risk.table.height = 0.5 ) #Useful when you have multiple groups 



surv_ext_test <- Surv(ext_test$DFS.time1, ext_test$DFS)
fit_ext_te <- survfit(surv_ext_test ~ ext_test$Risk, data= ext_test)

fit_ext_te_age <- survfit(surv_ext_test ~ (ext_test$age <median(ext_test$age)), data= ext_test)
fit_ext_te_gender <- survfit(surv_ext_test ~ ext_test$gender, data= ext_test)
fit_ext_te_stage <- survfit(surv_ext_test ~ ext_test$stage, data= ext_test)

ggsurvplot(fit_ext_te, pval = TRUE, censor = TRUE, conf.int = TRUE, palette = c("blue", "red"), break.time.by = 365, risk.table = TRUE, risk.table.height = 0.5 ) #Useful when you have multiple groups 

#ggsurvplot(fit_tr_age)
ggsurvplot(fit_ext_te_age, pval = TRUE, censor = TRUE,  break.time.by = 365, risk.table = TRUE, risk.table.height = 0.3 ) #Useful when you have multiple groups 

ggsurvplot(fit_ext_te_gender, pval = TRUE, censor = TRUE,  break.time.by = 365,  risk.table = TRUE, risk.table.height = 0.3 ) #Useful when you have multiple groups 

ggsurvplot(fit_ext_te_stage, pval = TRUE, censor = TRUE, conf.int = F,  break.time.by = 365, risk.table = TRUE, risk.table.height = 0.3 ) #Useful when you have multiple groups 

#pp_plot<- ggsurvplot(jj,data = data1,xlab = "Time in months", legend.labs = c(" >mean", " < mean"),legend=c(0.9,0.9), risk.table = TRUE,break.time.by = 20,risk.table.y.text.col = T, risk.table.y.text = FALSE, risk.table.height = 0.30)
pp_plot







