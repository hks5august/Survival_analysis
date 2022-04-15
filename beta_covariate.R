library(survival)
library(survminer)
library(dplyr)


############################ Load Training data #######################
train_mat<-read.csv("tcga_q_t",header =TRUE, sep = ",")
############################ Load validation data #######################
test1_mat<-read.csv("GSE14520_q_t",header =TRUE, sep = ",")
test2_mat<-read.csv("GSE76427_q_t",header =TRUE, sep = ",")

all_features_results<-read.csv("results/univariate_analysis/OS-final-HRmedian-quantright.csv",header =TRUE, sep = ",")

################################ Load selected feature files ##############################
features<- read.csv("results/OS/OS_Clin_both_CI_feature_names", header =TRUE, sep = ",", dec = ".",stringsAsFactors=FALSE)

features1<-features[,1]

################################ Make final files with the selected features (genes here) from the above block ###########################################################
train_feature_mat<- cbind(train_mat["OS"],train_mat["OS.time"],train_mat[,features1])### Three blocks when files are loaded
test1_feature_mat<- cbind(test1_mat["OS"],test1_mat["OS.time"],test1_mat[,features1])
test2_feature_mat<- cbind(test2_mat["OS"],test2_mat["OS.time"],test2_mat[,features1])

#train_feature_mat<- cbind(train_mat["gender"],train_mat["age"],train_mat["OS"],train_mat["OS.time"],train_mat[,features])### Three blocks when files are loaded
#test1_feature_mat<- cbind(test_mat1["gender"],test_mat1["age"],test_mat1["OS"],test_mat1["OS.time"],test_mat1[,features])
#test2_feature_mat<- cbind(test_mat2["gender"],test_mat2["age"],test_mat2["OS"],test_mat2["OS.time"],test_mat2[,features])

################## results of selected features ### extract rows from the dataframe ####################
sel_features_results<-all_features_results[row.names(all_features_results) %in% row.names(features),]

############ remove where OS.time=NA ############

train_feature_mat<-subset(train_feature_mat,OS.time!="nan")
test1_feature_mat<-subset(test1_feature_mat,OS.time!="nan")
test2_feature_mat<-subset(test2_feature_mat,OS.time!="nan")


Lt<- data.frame(sel_features_results["Gene"]);
Lt<-t(Lt)
LUSC_C=train_feature_mat

E= length(LUSC_C)

PI=0
for(i in seq(from=3, to=E ,by=1))
{
  PI= PI+(LUSC_C[,i]*sel_features_results[,2])
}

LUSC_C$PI<-PI
########################################

surv_object <- Surv(time = LUSC_C$OS.time, event = LUSC_C$OS)


# ###survival analysis: fits cox ph model to find HR for PI
fit <- survfit(surv_object~(LUSC_C$PI>mean(LUSC_C$PI)), data=LUSC_C)
summary(fit)
fit.coxph <- coxph(surv_object ~(LUSC_C$PI>mean(LUSC_C$PI)), data=LUSC_C)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
summary(fit.coxph)
fit$n[1]
fit$n[2]

################################## Without PI i.e. only actual values ############
PI_w=0

for(i in seq(from=3, to=E ,by=1))
{
  PI_w= PI_w+(LUSC_C[,i]*1)
}


LUSC_C$PIw<-PI_w

surv_object_w <- Surv(time = LUSC_C$OS.time, event = LUSC_C$OS)


# ###survival analysis: fits cox ph model to find HR for PI
fit_w <- survfit(surv_object_w~(LUSC_C$PIw>mean(LUSC_C$PIw)), data=LUSC_C)
summary(fit_w)
fit.coxph_w <- coxph(surv_object_w ~(LUSC_C$PIw>mean(LUSC_C$PIw)), data=LUSC_C)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
summary(fit.coxph_w)
fit$n[1]
fit$n[2]

