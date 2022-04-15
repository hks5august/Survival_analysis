library(survival)
library(survminer)
library(dplyr)


############################ Load Training data #######################
train_mat<-read.csv("../../../tcga_q_t",header =TRUE, sep = ",")
############################ Load validation data #######################
test1_mat<-read.csv("../../../GSE14520_q_t",header =TRUE, sep = ",")
test1_mat<-subset(test1_mat,OS.time!="nan")

test2_mat<-read.csv("../../../GSE76427_q_t",header =TRUE, sep = ",")
test2_mat<-subset(test2_mat,OS.time!="nan")


pred_test1<-read.csv("OS_time_Prediction_VH_vimp_valid1_lasso",sep = ",",header = TRUE)
pred_test2<-read.csv("OS_time_Prediction_VH_vimp_valid2_lasso",sep = ",",header = TRUE)
all_features_results<-read.csv("../../../results/univariate_analysis/OS-final-HRmedian-quantright.csv",header =TRUE, sep = ",")

################################ Load selected feature files ##############################
features<- read.csv("../../../results/OS/OS_Clin_both_MD_feature_names", header =TRUE, sep = ",", dec = ".",stringsAsFactors=FALSE)

features1<-features[,1]

################################ Make final files with the selected features (genes here) from the above block ###########################################################
train_feature_mat<- cbind(train_mat["OS"],train_mat["OS.time"],train_mat[,features1])### Three blocks when files are loaded
test1_feature_mat<- cbind(test1_mat["OS"],test1_mat["OS.time"],test1_mat[,features1])
test2_feature_mat<- cbind(test2_mat["OS"],test2_mat["OS.time"],test2_mat[,features1])

pred_t1<-cbind(test1_mat["OS"],test1_mat["OS.time"],pred_test1$OS_time_valid1_lasso)
pred_t2<-cbind(test2_mat["OS"],test2_mat["OS.time"],pred_test2$OS_time_valid2_lasso)



#train_feature_mat<- cbind(train_mat["gender"],train_mat["age"],train_mat["OS"],train_mat["OS.time"],train_mat[,features])### Three blocks when files are loaded
#test1_feature_mat<- cbind(test_mat1["gender"],test_mat1["age"],test_mat1["OS"],test_mat1["OS.time"],test_mat1[,features])
#test2_feature_mat<- cbind(test_mat2["gender"],test_mat2["age"],test_mat2["OS"],test_mat2["OS.time"],test_mat2[,features])

################## results of selected features ### extract rows from the dataframe ####################
sel_features_results<-all_features_results[row.names(all_features_results) %in% row.names(features),]

############ remove where OS.time=NA ############

train_feature_mat<-subset(train_feature_mat,OS.time!="nan")
test1_feature_mat<-subset(test1_feature_mat,OS.time!="nan")
test2_feature_mat<-subset(test2_feature_mat,OS.time!="nan")

##################################### Training dataset1 ##########################################
############################ Survival analysis using PI i.e. of selected genes features ###########

#Lt<- data.frame(sel_features_results["Gene"]);
#Lt<-t(Lt)
LUSC_C=train_feature_mat
#LUSC_C=test1_feature_mat
#LUSC_C=test2_feature_mat


E= length(LUSC_C)
PI=0
for(i in seq(from=3, to=E ,by=1))
{
  PI= PI+(LUSC_C[,i]*sel_features_results[i-2,2])
}

LUSC_C$PI<-PI
########################################
surv_object <- Surv(time = LUSC_C$OS.time, event = LUSC_C$OS)
# ###survival analysis: fits cox ph model to find HR for PI
fit <- survfit(surv_object~(LUSC_C$PI>median(LUSC_C$PI)), data=LUSC_C)
summary(fit)

fit.coxph <- coxph(surv_object ~(LUSC_C$PI>median(LUSC_C$PI)), data=LUSC_C)
summary(fit.coxph)
#ggsurvplot(fit)
summary(fit.coxph)
fit$n[1]
fit$n[2]

first <- coef(summary(fit.coxph))
CI<-exp(confint(fit.coxph))
#as.matrix(first)

#check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
#if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
#{

write.table(cbind("Beta","HR","CI-lower","CI-upper","P-value","GP1","GP2","Hr-Inv-lst","CI-lower","CI-upper","Concordance","Std_Error"),
            file="final_results_only_sel_features_with_PI_training.csv",row.names=F,col.names=F,sep = ',');
write.table(cbind(first[1],first[2],CI[1],CI[2],first[5],fit$n[1],fit$n[2],1/first[2],1/CI[1],1/CI[2],fit.coxph$concordance[6],fit.coxph$concordance[7]),
            file="final_results_only_sel_features_with_PI_training.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
#}


############################ Survival analysis using selected genes features ###########

################################## Without PI i.e. only actual values ############
PI_w=0

for(i in seq(from=3, to=E ,by=1))
{
  PI_w= PI_w+(LUSC_C[,i]*1)
}


LUSC_C$PIw<-PI_w

surv_object_w <- Surv(time = LUSC_C$OS.time, event = LUSC_C$OS)


# ###survival analysis: fits cox ph model to find HR for PI
#fit_w <- survfit(surv_object_w~(LUSC_C$PIw>mean(LUSC_C$PIw)), data=LUSC_C)
#summary(fit_w)
#fit.coxph_w <- coxph(surv_object_w ~(LUSC_C$PIw>mean(LUSC_C$PIw)), data=LUSC_C)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
#summary(fit.coxph_w)

fit_w1 <- survfit(surv_object_w~(LUSC_C$PIw>median(LUSC_C$PIw)), data=LUSC_C)
summary(fit_w1)
fit.coxph_w1 <- coxph(surv_object_w ~(LUSC_C$PIw>median(LUSC_C$PIw)), data=LUSC_C)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
summary(fit.coxph_w1)
first1 <- coef(summary(fit.coxph_w1))
CI1<-exp(confint(fit.coxph_w1))
#as.matrix(first)
fit1$n[1]
fit1$n[2]
#check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
#if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
#{
write.table(cbind("Beta","HR","CI-lower","CI-upper","P-value","GP1","GP2","Hr-Inv-lst","CI-lower","CI-upper","Concordance","Std_Error"),
            file="final_results_only_sel_features_without_PI_training.csv",row.names=F,col.names=F,sep = ',');
write.table(cbind(first1[1],first1[2],CI1[1],CI1[2],first1[5],fit_w1$n[1],fit_w1$n[2],1/first1[2],1/CI1[1],1/CI1[2],fit.coxph_w1$concordance[6],fit.coxph_w1$concordance[7]),
            file="final_results_only_sel_features_without_PI_training.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
#}


##################################################################################################


############################ Survival analysis using PI i.e. of selected genes features ###########
############################ Validation dataset1 #################################################

#Lt<- data.frame(sel_features_results["Gene"]);
#Lt<-t(Lt)
#LUSC_C=train_feature_mat
LUSC_C=test1_feature_mat
#LUSC_C=test2_feature_mat


E= length(LUSC_C)
PI=0
for(i in seq(from=3, to=E ,by=1))
{
  PI= PI+(LUSC_C[,i]*sel_features_results[i-2,2])
}

LUSC_C$PI<-PI
########################################
surv_object <- Surv(time = LUSC_C$OS.time, event = LUSC_C$OS)
# ###survival analysis: fits cox ph model to find HR for PI
fit <- survfit(surv_object~(LUSC_C$PI>median(LUSC_C$PI)), data=LUSC_C)
summary(fit)

fit.coxph <- coxph(surv_object ~(LUSC_C$PI>median(LUSC_C$PI)), data=LUSC_C)
summary(fit.coxph)
#ggsurvplot(fit)
summary(fit.coxph)
fit$n[1]
fit$n[2]

first <- coef(summary(fit.coxph))
CI<-exp(confint(fit.coxph))
#as.matrix(first)

#check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
#if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
#{

write.table(cbind("Beta","HR","CI-lower","CI-upper","P-value","GP1","GP2","Hr-Inv-lst","CI-lower","CI-upper","Concordance","Std_Error"),
            file="final_results_only_sel_features_with_PI_valid1.csv",row.names=F,col.names=F,sep = ',');
write.table(cbind(first[1],first[2],CI[1],CI[2],first[5],fit$n[1],fit$n[2],1/first[2],1/CI[1],1/CI[2],fit.coxph$concordance[6],fit.coxph$concordance[7]),
            file="final_results_only_sel_features_with_PI_valid1.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
#}


############################ Survival analysis using selected genes features ###########

################################## Without PI i.e. only actual values ############
PI_w=0

for(i in seq(from=3, to=E ,by=1))
{
  PI_w= PI_w+(LUSC_C[,i]*1)
}


LUSC_C$PIw<-PI_w

surv_object_w <- Surv(time = LUSC_C$OS.time, event = LUSC_C$OS)


# ###survival analysis: fits cox ph model to find HR for PI
#fit_w <- survfit(surv_object_w~(LUSC_C$PIw>mean(LUSC_C$PIw)), data=LUSC_C)
#summary(fit_w)
#fit.coxph_w <- coxph(surv_object_w ~(LUSC_C$PIw>mean(LUSC_C$PIw)), data=LUSC_C)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
#summary(fit.coxph_w)

fit_w1 <- survfit(surv_object_w~(LUSC_C$PIw>median(LUSC_C$PIw)), data=LUSC_C)
summary(fit_w1)
fit.coxph_w1 <- coxph(surv_object_w ~(LUSC_C$PIw>median(LUSC_C$PIw)), data=LUSC_C)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
summary(fit.coxph_w1)
first1 <- coef(summary(fit.coxph_w1))
CI1<-exp(confint(fit.coxph_w1))
#as.matrix(first)
fit1$n[1]
fit1$n[2]
#check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
#if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
#{
write.table(cbind("Beta","HR","CI-lower","CI-upper","P-value","GP1","GP2","Hr-Inv-lst","CI-lower","CI-upper","Concordance","Std_Error"),
            file="final_results_only_sel_features_without_PI_valid1.csv",row.names=F,col.names=F,sep = ',');
write.table(cbind(first1[1],first1[2],CI1[1],CI1[2],first1[5],fit_w1$n[1],fit_w1$n[2],1/first1[2],1/CI1[1],1/CI1[2],fit.coxph_w1$concordance[6],fit.coxph_w1$concordance[7]),
            file="final_results_only_sel_features_without_PI_valid1.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
#}


##################################################################################################
######################################### on validation dataset2 #################################

############################ Survival analysis using PI i.e. of selected genes features ###########

#Lt<- data.frame(sel_features_results["Gene"]);
#Lt<-t(Lt)
#LUSC_C=train_feature_mat
#LUSC_C=test1_feature_mat
LUSC_C=test2_feature_mat


E= length(LUSC_C)
PI=0
for(i in seq(from=3, to=E ,by=1))
{
  PI= PI+(LUSC_C[,i]*sel_features_results[i-2,2])
}

LUSC_C$PI<-PI
########################################
surv_object <- Surv(time = LUSC_C$OS.time, event = LUSC_C$OS)
# ###survival analysis: fits cox ph model to find HR for PI
fit <- survfit(surv_object~(LUSC_C$PI>median(LUSC_C$PI)), data=LUSC_C)
summary(fit)

fit.coxph <- coxph(surv_object ~(LUSC_C$PI>median(LUSC_C$PI)), data=LUSC_C)
summary(fit.coxph)
#ggsurvplot(fit)
summary(fit.coxph)
fit$n[1]
fit$n[2]

first <- coef(summary(fit.coxph))
CI<-exp(confint(fit.coxph))
#as.matrix(first)

#check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
#if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
#{

write.table(cbind("Beta","HR","CI-lower","CI-upper","P-value","GP1","GP2","Hr-Inv-lst","CI-lower","CI-upper","Concordance","Std_Error"),
            file="final_results_only_sel_features_with_PI_valid2.csv",row.names=F,col.names=F,sep = ',');
write.table(cbind(first[1],first[2],CI[1],CI[2],first[5],fit$n[1],fit$n[2],1/first[2],1/CI[1],1/CI[2],fit.coxph$concordance[6],fit.coxph$concordance[7]),
            file="final_results_only_sel_features_with_PI_valid2.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
#}


############################ Survival analysis using selected genes features ###########

################################## Without PI i.e. only actual values ############
PI_w=0

for(i in seq(from=3, to=E ,by=1))
{
  PI_w= PI_w+(LUSC_C[,i]*1)
}


LUSC_C$PIw<-PI_w

surv_object_w <- Surv(time = LUSC_C$OS.time, event = LUSC_C$OS)


# ###survival analysis: fits cox ph model to find HR for PI
#fit_w <- survfit(surv_object_w~(LUSC_C$PIw>mean(LUSC_C$PIw)), data=LUSC_C)
#summary(fit_w)
#fit.coxph_w <- coxph(surv_object_w ~(LUSC_C$PIw>mean(LUSC_C$PIw)), data=LUSC_C)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
#summary(fit.coxph_w)

fit_w1 <- survfit(surv_object_w~(LUSC_C$PIw>median(LUSC_C$PIw)), data=LUSC_C)
summary(fit_w1)
fit.coxph_w1 <- coxph(surv_object_w ~(LUSC_C$PIw>median(LUSC_C$PIw)), data=LUSC_C)
#ggsurvplot(fit)
# write.csv(LUSC_C, file = "/Users/sherrybhalla/Desktop/LUNG-Cancer/LUNG-Apop-HRmean-PImean.csv",row.names = FALSE)
summary(fit.coxph_w1)
first1 <- coef(summary(fit.coxph_w1))
CI1<-exp(confint(fit.coxph_w1))
#as.matrix(first)
fit1$n[1]
fit1$n[2]
#check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
#if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
#{
write.table(cbind("Beta","HR","CI-lower","CI-upper","P-value","GP1","GP2","Hr-Inv-lst","CI-lower","CI-upper","Concordance","Std_Error"),
            file="final_results_only_sel_features_without_PI_valid2.csv",row.names=F,col.names=F,sep = ',');
write.table(cbind(first1[1],first1[2],CI1[1],CI1[2],first1[5],fit_w1$n[1],fit_w1$n[2],1/first1[2],1/CI1[1],1/CI1[2],fit.coxph_w1$concordance[6],fit.coxph_w1$concordance[7]),
            file="final_results_only_sel_features_without_PI_valid2.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
#}





#################################################################################################



############################ Survival analysis using predicted OS time as features ###########

surv_object2 <- Surv(time = pred_t1$OS.time, event = pred_t1$OS)



#fit2 <- survfit(surv_object2~(pred_t1$`pred_test1$OS_time_valid1_rf`>mean(pred_t1$`pred_test1$OS_time_valid1_rf`)), data=LUSC_C)
#summary(fit2)
fit2 <- survfit(surv_object2~(pred_t1$`pred_test1$OS_time_valid1_lasso`>median(pred_t1$`pred_test1$OS_time_valid1_lasso`)), data=LUSC_C)
summary(fit2)

#ggsurvplot(fit)
#ggsurvplot(fit1)

#fit.coxph2 <- coxph(surv_object2 ~(pred_t1$`pred_test1$OS_time_valid1_rf`>mean(pred_t1$`pred_test1$OS_time_valid1_rf`)), data=LUSC_C)
#summary(fit.coxph2)
fit.coxph2 <- coxph(surv_object2 ~(pred_t1$`pred_test1$OS_time_valid1_lasso`>median(pred_t1$`pred_test1$OS_time_valid1_lasso`)), data=LUSC_C)
summary(fit.coxph2)

first2 <- coef(summary(fit.coxph2))
CI2<-exp(confint(fit.coxph2))
#as.matrix(first)
fit2$n[1]
fit2$n[2]
#check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
#if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
#{
write.table(cbind("Beta","HR","CI-lower","CI-upper","P-value","GP1","GP2","Hr-Inv-lst","CI-lower","CI-upper","Concordance","Std_Error"),
            file="final_results_with_predicted_OS_time_valid1.csv",row.names=F,col.names=F,sep = ',');
write.table(cbind(first2[1],first2[2],CI2[1],CI2[2],first2[5],fit2$n[1],fit2$n[2],1/first2[2],1/CI2[1],1/CI2[2],fit.coxph2$concordance[6],fit.coxph2$concordance[7]),
            file="final_results_with_predicted_OS_time_valid1.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
#}


############################ Survival analysis using predicted OS time as features of validation data 2 ###########

surv_object2 <- Surv(time = pred_t2$OS.time, event = pred_t2$OS)



#fit2 <- survfit(surv_object2~(pred_t1$`pred_test1$OS_time_valid1_rf`>mean(pred_t1$`pred_test1$OS_time_valid1_rf`)), data=LUSC_C)
#summary(fit2)
fit2 <- survfit(surv_object2~(pred_t2$`pred_test2$OS_time_valid2_lasso`>median(pred_t2$`pred_test2$OS_time_valid2_lasso`)), data=LUSC_C)
summary(fit2)

#ggsurvplot(fit)
#ggsurvplot(fit1)

#fit.coxph2 <- coxph(surv_object2 ~(pred_t1$`pred_test1$OS_time_valid1_rf`>mean(pred_t1$`pred_test1$OS_time_valid1_rf`)), data=LUSC_C)
#summary(fit.coxph2)
fit.coxph2 <- coxph(surv_object2 ~(pred_t2$`pred_test2$OS_time_valid2_lasso`>median(pred_t2$`pred_test2$OS_time_valid2_lasso`)), data=LUSC_C)
summary(fit.coxph2)

first2 <- coef(summary(fit.coxph2))
CI2<-exp(confint(fit.coxph2))
#as.matrix(first)
fit2$n[1]
fit2$n[2]
#check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
#if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
#{
write.table(cbind("Beta","HR","CI-lower","CI-upper","P-value","GP1","GP2","Hr-Inv-lst","CI-lower","CI-upper","Concordance","Std_Error"),
            file="final_results_with_predicted_OS_time_valid2.csv",row.names=F,col.names=F,sep = ',');
write.table(cbind(first2[1],first2[2],CI2[1],CI2[2],first2[5],fit2$n[1],fit2$n[2],1/first2[2],1/CI2[1],1/CI2[2],fit.coxph2$concordance[6],fit.coxph2$concordance[7]),
            file="final_results_with_predicted_OS_time_valid2.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
#}




