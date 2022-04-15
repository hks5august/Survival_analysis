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


ext_test <- read.csv(file="final_ext_te1.csv", header =TRUE, sep = ",", dec = ".", );
################## remove samples where OS.time or DFI time is NA ##########

ext_test<-subset(ext_test,DFS.time1!="nan")


surv_train <- Surv(train$DFS.time, train$DFS)
#fit_tr <- coxph(surv_train ~ group, data= train)
#fit_tr <- coxph(surv_train ~ train$pred, data= train)

fit_tr <- survfit(surv_train ~ train$Risk, data= train)
fit_tr_age <- survfit(surv_train ~ (train$age <median(train$age)), data= train)
fit_tr_gender <- survfit(surv_train ~ train$gender, data= train)
fit_tr_stage <- survfit(surv_train ~ train$stage, data= train)

ggsurvplot(fit_tr_stage, pval = TRUE, censor = TRUE, conf.int = F,  break.time.by = 365, xlab = "Time in days", risk.table = TRUE, risk.table.height = 0.3 ) #Useful when you have m

ggsurvplot(fit_tr_stage, pval = TRUE, censor = TRUE, conf.int = F, xlab = "Time in days", break.time.by = 365,  legend=c(0.8,0.8)) #Useful when you have m



######## Test Data ###########
surv_test <- Surv(test$DFS.time, test$DFS)
fit_te <- survfit(surv_test ~ test$Risk, data= test)

fit_te_age <- survfit(surv_test ~ (test$age <median(test$age)), data= test)
fit_te_gender <- survfit(surv_test ~ test$gender, data= test)
fit_te_stage <- survfit(surv_test ~ test$stage, data= test)



ggsurvplot(fit_te_stage, pval = TRUE, censor = TRUE, conf.int = F, xlab = "Time in days", break.time.by = 365, risk.table = TRUE, risk.table.height = 0.28 ) #Useful when you have m

ggsurvplot(fit_te_stage,  censor = TRUE, pval = TRUE,  conf.int = F, xlab = "Time in days", break.time.by = 365, legend=c(0.8,0.8)) #Useful when you have m




#### external validation #####
surv_ext_test <- Surv(ext_test$DFS.time1, ext_test$DFS)
fit_ext_te <- survfit(surv_ext_test ~ ext_test$Risk, data= ext_test)

fit_ext_te_age <- survfit(surv_ext_test ~ (ext_test$age <median(ext_test$age)), data= ext_test)
fit_ext_te_gender <- survfit(surv_ext_test ~ ext_test$gender, data= ext_test)
fit_ext_te_stage <- survfit(surv_ext_test ~ ext_test$stage, data= ext_test)



ggsurvplot(fit_ext_te_stage, pval = TRUE, censor = TRUE, conf.int = F, xlab = "Time in days", break.time.by = 365, risk.table = TRUE, risk.table.height = 0.28 ) #Useful when you have m

ggsurvplot(fit_ext_te_stage, pval = TRUE, censor = TRUE, conf.int = F,  xlab = "Time in days", break.time.by = 365, legend=c(0.8,0.8), size=1.5, font.tickslab = c(12, "bold", "black"), font.y = c(14, "bold", "black"), font.x = c(14, "bold", "black"), font.legend = c(14, "bold", "black"), font.label = c(14, "bold", "black"), test.for.trend = TRUE) #Useful when you have m



####
ggsurv<- ggsurvplot(fit_tr_stage,  censor = TRUE,  conf.int = F, xlab = "Time in days", break.time.by = 365, legend=c(0.8,0.8), size=1.5, font.tickslab = c(12, "bold", "black"), font.y = c(14, "bold", "black"), font.x = c(14, "bold", "black"), font.legend = c(14, "bold", "black"))
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 250, y = 0.18, # x and y coordinates of the text
                    label = "HR = 1.73 \n            P-value = <0.0001", size = 5,  font.label = c(14, "bold", "black"))
ggsurv                              
                                     
####
ggsurv<- ggsurvplot(fit_te_stage, censor = TRUE,  conf.int = F, xlab = "Time in days", break.time.by = 365, legend=c(0.8,0.8), size=1.5, font.tickslab = c(12, "bold", "black"), font.y = c(14, "bold", "black"), font.x = c(14, "bold", "black"), font.legend = c(14, "bold", "black"))
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 250, y = 0.18, # x and y coordinates of the text
                    label = "HR = 1.41 \n       P-value = 0.14", size = 5,  font.label = c(14, "bold", "black"))
ggsurv                              

####

ggsurv<- ggsurvplot(fit_ext_te_stage, censor = TRUE,  conf.int = F, xlab = "Time in days", break.time.by = 365, legend=c(0.8,0.8), size=1.5, font.tickslab = c(12, "bold", "black"), font.y = c(14, "bold", "black"), font.x = c(14, "bold", "black"), font.legend = c(14, "bold", "black"))
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 250, y = 0.28, # x and y coordinates of the text
                    label = "HR = 1.73 \n            P-value = <0.0001", size = 5,  font.label = c(14, "bold", "black"))
ggsurv                              




################################
surv_train <- Surv(train$DFS.time, train$DFS)
fit_tr <- survfit(surv_train ~ train$Risk_prediction, data= train)

ggsurv<- ggsurvplot(fit_tr, censor = TRUE,  conf.int = F, xlab = "Time in days", break.time.by = 365, legend.labs = c(" High-Risk", "Low-Risk"), legend=c(0.8,0.8), size=1.5, font.tickslab = c(12, "bold", "black"), font.y = c(14, "bold", "black"), font.x = c(14, "bold", "black"), font.legend = c(14, "bold", "black"))
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 250, y = 0.28, # x and y coordinates of the text
                    label = "HR = 2.91 \n            P-value = <0.0001", size = 5,  font.label = c(14, "bold", "black"))
ggsurv                              

#################################

surv_test <- Surv(test$DFS.time, test$DFS)
fit_te <- survfit(surv_test ~ test$Risk, data= test)

ggsurv<- ggsurvplot(fit_te, censor = TRUE,  conf.int = F, xlab = "Time in days", break.time.by = 365, legend.labs = c(" High-Risk", "Low-Risk"),legend=c(0.8,0.8), size=1.5, font.tickslab = c(12, "bold", "black"), font.y = c(14, "bold", "black"), font.x = c(14, "bold", "black"), font.legend = c(14, "bold", "black"))
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 250, y = 0.28, # x and y coordinates of the text
                    label = "HR = 3.65 \n           P-value = 0.0014", size = 5,  font.label = c(14, "bold", "black"))
ggsurv                              


#################################
surv_ext_test <- Surv(ext_test$DFS.time1, ext_test$DFS)
fit_ext_te <- survfit(surv_ext_test ~ ext_test$Risk_prediction, data= ext_test)


ggsurv<- ggsurvplot(fit_ext_te, censor = TRUE,  conf.int = F, xlab = "Time in days", break.time.by = 365,legend.labs = c(" High-Risk", "Low-Risk"), legend=c(0.8,0.8), size=1.5, font.tickslab = c(12, "bold", "black"), font.y = c(14, "bold", "black"), font.x = c(14, "bold", "black"), font.legend = c(14, "bold", "black"))
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 250, y = 0.28, # x and y coordinates of the text
                    label = "HR = 1.79 \n         P-value = 0.0025", size = 5,  font.label = c(14, "bold", "black")  )
ggsurv                              



