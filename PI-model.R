library(survival)
library(survminer)
library(dplyr)
library(survcomp)
library("survminer")
library(survAUC)


set.seed(7)

train_mat<-read.csv("4features_train_mat",header =TRUE, sep = ",")
############################ Load validation data #######################
test_mat<-read.csv("4features_test_mat",header =TRUE, sep = ",")


############ remove where OS.time=NA ############

train<-subset(train_mat,DFS.time!="nan")
test<-subset(test_mat,DFS.time!="nan")

train<-subset(train,DFS.time!="nan")
test<-subset(test,DFS.time!="nan")

LUSC_C=train

PI=0

PI= PI+LUSC_C$hsa.mir.3936*0.48+LUSC_C$SUZ12*0.39+LUSC_C$cg18465072*0.45+LUSC_C$cg22852503*(-0.45)
#PI= PI+LUSC_C$SUZ12*0.39
PI
LUSC_C$PI<-PI

surv_object <- Surv(time = LUSC_C$DFS.time, event = LUSC_C$DFS)


# ###survival analysis: fits cox ph model to find HR for PI
#fit <- survfit(surv_object~(LUSC_C$PI>mean(LUSC_C$PI)), data=LUSC_C)
#summary(fit)
################ KM Plots ################
##fit_tr <- survfit(surv_train ~ ((LUSC_C$PI>mean(LUSC_C$PI)), data= train)
#ggsurvplot(fit, pval = TRUE, censor = TRUE,  break.time.by = 365, risk.table = TRUE, risk.table.height = 0.25 ) #Useful when you have multiple
#ggsurvplot(fit, pval = TRUE, censor = TRUE, conf.int = F,  break.time.by = 365, xlab = "Time in days", risk.table = TRUE, risk.table.height = 0.3 )

fit <- survfit(surv_object~(LUSC_C$PI<mean(LUSC_C$PI)), data=LUSC_C)
#ggsurvplot(fit, pval = TRUE, censor = TRUE, conf.int = F, xlab = "Time in days", break.time.by = 365,  legend=c(0.8,0.8)) #Useful when you have m

ggsurv<- ggsurvplot(fit, censor = TRUE,  conf.int = F, xlab = "Time in days", break.time.by = 365, legend.labs = c(" High-Risk or PI > mean(PI)", "Low-Risk or PI < mean(PI)"), legend=c(0.8,0.8), size=1.5, font.tickslab = c(12, "bold", "black"), font.y = c(14, "bold", "black"), font.x = c(14, "bold", "black"), font.legend = c(14, "bold", "black"))
ggsurv$plot <- ggsurv$plot+
  ggplot2::annotate("text",
           x = 250, y = 0.28, # x and y coordinates of the text
                    label = "HR = 2.37 \n            P-value = <0.0001", size = 5,  font.label = c(14, "bold", "black"))
ggsurv



##### COXPH model #####


fit.coxph <- coxph(surv_object ~(LUSC_C$PI>mean(LUSC_C$PI)), data=LUSC_C)
#fit.coxph <- coxph(surv_object ~(LUSC_C$PI>median(LUSC_C$PI)), data=LUSC_C)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
summary(fit.coxph)
fit$n[1]
fit$n[2]


#tr_CI<-concordance.index(x=train$Ridge.Classifier, surv.time=train$DFS.time, surv.event=train$DFS, method="noether")

tr_CI<-concordance.index(x=(LUSC_C$PI>mean(LUSC_C$PI)), surv.time=train$DFS.time, surv.event=train$DFS, method="noether")

tr_CI


HR_tr<-hazard.ratio(x=(LUSC_C$PI>mean(LUSC_C$PI)), surv.time=train$DFS.time, surv.event=train$DFS, alpha = 0.05, method.test = c("logrank"), na.rm = T)

HR_tr


######## Brier Score ####

dd_tr <- data.frame("time"=train$DFS.time, "event"=train$DFS, "score"=(LUSC_C$PI>mean(LUSC_C$PI)))


tr_Brier_score <- sbrier.score2proba(data.tr=dd_tr, data.ts=dd_tr, method="cox")

tr_Brier_score$bsc.integrated


###Without coefficients

PI_w=0

PI_w= PI_w+LUSC_C$hsa.mir.3936+LUSC_C$SUZ12+LUSC_C$cg18465072+LUSC_C$cg22852503
PI_w
LUSC_C$PI_w<-PI_w

surv_object <- Surv(time = LUSC_C$DFS.time, event = LUSC_C$DFS)


# ###survival analysis: fits cox ph model to find HR for PI
fit1 <- survfit(surv_object~(LUSC_C$PI_w>mean(LUSC_C$PI_w)), data=LUSC_C)
summary(fit1)
fit1.coxph <- coxph(surv_object ~(LUSC_C$PI_w>mean(LUSC_C$PI_w)), data=LUSC_C)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
summary(fit1.coxph)
fit1$n[1]
fit1$n[2]



########### test Data #####
LUSC_t=test

PI_t=0

PI_t= PI_t+LUSC_t$hsa.mir.3936*0.48+LUSC_t$SUZ12*0.39+LUSC_t$cg18465072*0.45+LUSC_t$cg22852503*(-0.45)
#PI_t= PI_t+LUSC_t$SUZ12*0.39
PI_t
LUSC_t$PI<-PI_t

surv_object_t <- Surv(time = LUSC_t$DFS.time, event = LUSC_t$DFS)


# ###survival analysis: fits cox ph model to find HR for PI
fit_t <- survfit(surv_object_t~(LUSC_t$PI<mean(LUSC_t$PI)), data=LUSC_t)
summary(fit_t)


##### KM plot ####

#ggsurv<- ggsurvplot(fit_t, censor = TRUE,  conf.int = F, xlab = "Time in days", break.time.by = 365, legend.labs = c(" High-Risk", "Low-Risk"), legend=c(0.8,0.8), size=1.5, font.tickslab = c(12, "bold", "black"), font.y = c(14, "bold", "black"), font.x = c(14, "bold", "black"), font.legend = c(14, "bold", "black"))


#ggsurvplot(fit_t, pval = TRUE, censor = TRUE, conf.int = F,  break.time.by = 365, xlab = "Time in days", risk.table = TRUE, risk.table.height = 0.3 )

#ggsurvplot(fit_t, pval = TRUE, censor = TRUE, conf.int = F, xlab = "Time in days", break.time.by = 365,  legend=c(0.8,0.8)) #Useful when you have m

ggsurv<- ggsurvplot(fit_t, censor = TRUE,  conf.int = F, xlab = "Time in days", break.time.by = 365, legend.labs = c(" High-Risk or PI > mean(PI)", "Low-Risk or PI < mean(PI)"), legend=c(0.8,0.8), size=1.5, font.tickslab = c(12, "bold", "black"), font.y = c(14, "bold", "black"), font.x = c(14, "bold", "black"), font.legend = c(14, "bold", "black"))
ggsurv$plot <- ggsurv$plot+
  ggplot2::annotate("text",
                    x = 250, y = 0.28, # x and y coordinates of the text
                    label = "HR = 2.54 \n         P-value = 0.015", size = 5,  font.label = c(14, "bold", "black"))
ggsurv


##### COCPH model ####

fit_t.coxph <- coxph(surv_object_t ~(LUSC_t$PI>mean(LUSC_t$PI)), data=LUSC_t)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
summary(fit_t.coxph)
fit_t$n[1]
fit_t$n[2]

te_CI<-concordance.index(x=(LUSC_t$PI>mean(LUSC_t$PI)), surv.time=test$DFS.time, surv.event=test$DFS, method="noether")

te_CI

HR_te<-hazard.ratio(x=(LUSC_t$PI>mean(LUSC_t$PI)), surv.time=test$DFS.time, surv.event=test$DFS, alpha = 0.05, method.test = c("logrank"), na.rm = T)

HR_te


######## Brier Score ####

dd_te <- data.frame("time"=test$DFS.time, "event"=test$DFS, "score"=(LUSC_t$PI>mean(LUSC_t$PI)))


te_Brier_score <- sbrier.score2proba(data.tr=dd_te, data.ts=dd_te, method="cox")

te_Brier_score$bsc.integrated
############## 


###Without coefficients

PI_w_t=0

PI_w_t= PI_w_t+LUSC_t$hsa.mir.3936+LUSC_t$SUZ12+LUSC_t$cg18465072+LUSC_t$cg22852503
PI_w_t
LUSC_t$PI_w_t<-PI_w_t

surv_object <- Surv(time = LUSC_t$DFS.time, event = LUSC_t$DFS)


# ###survival analysis: fits cox ph model to find HR for PI
fit1_t <- survfit(surv_object~(LUSC_t$PI_w_t>mean(LUSC_t$PI_w_t)), data=LUSC_t)
summary(fit1)
fit1.coxph_t <- coxph(surv_object ~(LUSC_t$PI_w_t>mean(LUSC_t$PI_w_t)), data=LUSC_t)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
summary(fit1.coxph_t)
fit1$n[1]
fit1$n[2]



##############




