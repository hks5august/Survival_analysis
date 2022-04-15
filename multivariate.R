library(survival)
library(survminer)
library(dplyr)

############################ Load Training data #######################
mat<-read.csv("GSE3",header =TRUE, sep = "\t")


############################## create survival object for OS #############
surv_object_OS <- Surv(time = mat$Survival.months, event = mat$Survival.status)


############################## create survival object for RFS #############
surv_object_RFS <- Surv(time = mat$Recurr.months, event = mat$Recurr.status)


############################## Extract features on which multivariate performed ########

Age<-mat$Age
Gender<-mat$Gender
Stage<-mat$TNM.staging
FCN3<-mat$FCN3
CLEC1B<-mat$CLEC1B
PRC1<-mat$PRC1
comb<-mat$PRC1*0.428+mat$FCN3*(-0.463)+mat$CLEC1B*(-0.595)

############################## multivariate analysis #########################
cox_OS<- coxph(surv_object_OS ~Age+Gender+Stage+FCN3+CLEC1B+PRC1,data=mat)

####### plot ##########
ggforest(cox_OS,data=mat)
############### multivariate analysis on RFS ############
cox_RFS<- coxph(surv_object_RFS ~Age+Gender+Stage+FCN3+CLEC1B+PRC1,data=mat)

####### plot ##########
ggforest(cox_RFS,data=mat)

