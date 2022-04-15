library(survival)
library(survminer)
library(dplyr)

#### input data ####### features should be in columns and samples should be in Rows,  ###

input_data <- read.csv(file='tcga_q_t.csv',header =TRUE, sep = ",", dec = ".");
################## remove samples where OS.time or DFI time is NA ##########
input_data1<-subset(input_data,OS.time!="nan")

input_data2<-subset(input_data,DFI.time!="nan")
############# For overall Survival ##############


write.table(cbind("Gene","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="results/OS-final-HRmedian-quantright.csv",row.names=F,col.names=F,sep = ',');

################ Here features to compute survival start from 6th column onwards ##############

for(i in seq(from=6, to=length(input_data1), by=1))
{
  surv_object <- Surv(time = input_data1$OS.time, event = input_data$OS)
  
  #survival analysis: fits cox ph model to find HR for mean cut
  fit1 <- survfit(surv_object~(input_data1[,i])>(median(input_data[1,i])), data=input_data1);
  summary(fit1);
  fit1.coxph <- coxph(surv_object ~ (input_data1[,i])>(median(input_data1[,i])), data = input_data1)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  #as.matrix(first)
  
  #check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
  if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
  {write.table(cbind(colnames(input_data1[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="results/OS-final-HRmedian-quantright.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
  }
}


############# For Disease Free Survival ##############


write.table(cbind("Gene","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="results/DFS-final-HRmedian-quantright.csv",row.names=F,col.names=F,sep = ',');
################ Here features to compute survival start from 6th column onwards ##############

for(i in seq(from=6, to=length(input_data2), by=1))
{
  surv_object <- Surv(time = input_data2$DFI.time, event = input_data2$DFI)
  
  #survival analysis: fits cox ph model to find HR for mean cut
  fit1 <- survfit(surv_object~(input_data2[,i])>(median(input_data2[,i])), data=input_data2);
  summary(fit1);
  fit1.coxph <- coxph(surv_object ~ (input_data2[,i])>(median(input_data2[,i])), data = input_data2)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  #as.matrix(first)
  
  #check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
  if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
  {write.table(cbind(colnames(input_data2[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="results/DFS-final-HRmedian-quantright.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
  }
}

######### on the basis of mean cut off ##########

############# For overall Survival ##############


write.table(cbind("Gene","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="results/OS-final-HRmean-quantright.csv",row.names=F,col.names=F,sep = ',');

################ Here features to compute survival start from 6th column onwards ##############
for(i in seq(from=6, to=length(input_data1), by=1))
{
  surv_object <- Surv(time = input_data1$OS.time, event = input_data1$OS)
  
  #survival analysis: fits cox ph model to find HR for mean cut
  fit1 <- survfit(surv_object~(input_data1[,i])>(mean(input_data1[,i])), data=input_data1);
  summary(fit1);
  fit1.coxph <- coxph(surv_object ~ (input_data1[,i])>(mean(input_data1[,i])), data = input_data1)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  #as.matrix(first)
  
  #check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
  if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
  {write.table(cbind(colnames(input_data1[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="results/OS-final-HRmean-quantright.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
  }
}






############# For Disease Free Survival ##############


write.table(cbind("Gene","Beta","HR","P-value","GP1","GP2","Hr-Inv-lst","Concordance","Std_Error"),
            file="results/DFS-final-HRmean-quantright.csv",row.names=F,col.names=F,sep = ',');

################ Here features to compute survival start from 6th column onwards ##############
for(i in seq(from=6, to=length(input_data2), by=1))
{
  surv_object <- Surv(time = input_data2$DFI.time, event = input_data2$DFI)
  
  #survival analysis: fits cox ph model to find HR for mean cut
  fit1 <- survfit(surv_object~(input_data2[,i])>(mean(input_data2[,i])), data=input_data2);
  summary(fit1);
  fit1.coxph <- coxph(surv_object ~ (input_data2[,i])>(mean(input_data2[,i])), data = input_data2)
  # summary(fit1.coxph);
  first <- coef(summary(fit1.coxph))
  #as.matrix(first)
  
  #check whether the pvalue is significant and HR is more than 20 (only bclxl had HR=20.8)
  if((first[5]<=0.05)&&(!is.na(first[5]))&&(!is.na(first[2])))
  {write.table(cbind(colnames(input_data2[i]),first[1],first[2],first[5],fit1$n[1],fit1$n[2],1/first[2],fit1.coxph$concordance[6],fit1.coxph$concordance[7]),
               file="results/DFS-final-HRmean-quantright.csv",row.names=F,col.names=F,sep = ',',append = T);#output file
  }
}

