### load libraries 
library("rpart")
library("ipred")
library("survival")
library("survminer")
library("DStree")
library("party")
library("plyr")
library("RANN")
library("dplyr")
library("RSNNS")
library("kernlab")
# library("neuralnet")
# library("xgboost")
# library("brnn")
library("randomForestSRC")
library("randomForest")
library("mlbench") ### for MAE calculation
library("keras")
library("glmnet")
library("caret")


############ load complete data file ############
input_data <- read.csv(file='tcga_q',header =TRUE, sep = ",", dec = ".");

############# load significant HR results file obtained from univariate analysis ##############

result_dfs<-read.csv(file='results/univariate_analysis/DFS-final-HRmedian-quantright.csv',header =TRUE, sep = ",", stringsAsFactors = FALSE);

############## Extract survival favorable and survival unfav genes #############

Surv_unfav <-subset(result_dfs, HR>1.2)                  # High gene expression is responsible for low survival rate [Survival_unfavourable].
Surv_fav <-subset(result_dfs, Hr.Inv.lst>1.2)          # High gene expression is responsible for high survival rate [Survival_favourable].

Surv_unfav<- Surv_unfav[,"Gene"]
Surv_unfav_g <- names(input_data)[(names(input_data) %in% Surv_unfav)]

write.table(Surv_unfav_g,file="results/univariate_analysis/list_DFS_surv_unfav_genes",row.names=F,col.names=F,sep = ',')
DFS_unfav<- input_data[,Surv_unfav_g]


Surv_fav<- Surv_fav[,"Gene"]
Surv_fav_g <- names(input_data)[(names(input_data) %in% Surv_fav)]

write.table(Surv_fav_g,file="results/univariate_analysis/list_DFS_surv_fav_genes",row.names=F,col.names=F,sep = ',')
DFS_fav<- input_data[,Surv_fav_g]
 ###### combination of fav and unfav ########
DFS_comb<- cbind(DFS_unfav,DFS_fav)    

DFS_Clin_unfav<- cbind(input_data["DFI"],input_data["DFI.time"],DFS_unfav)          # Unfavourable survival time gene expressions dataframe with clinical information
DFS_Clin_fav<- cbind(input_data["DFI"],input_data["DFI.time"],DFS_fav)            # Favourable survival time gene expressions dataframe with clinical information
DFS_Clin_both_fav_unfav<- cbind(input_data["DFI"],input_data["DFI.time"],DFS_comb)           # Favourable and unfavourable survival time gene expressions dataframe with clinical information


set.seed(17)

####### Feature selection form the data files... Below three nh code snippets take half an hour each to run on this dataset

######################## feature selection for favorable features ####################################

##### variable selection with "VH" or variable hunting method #####
DFS_Clin_fav_VH_features<- var.select.rfsrc(Surv(DFI.time, DFI) ~ ., DFS_Clin_fav,method = "vh", nrep = 100, nstep = 5)

#### Extract names of variables #####

DFS_Clin_fav_VH_feature_names <- DFS_Clin_fav_VH_features$topvars

write.table(DFS_Clin_fav_VH_feature_names,file="DFS_Clin_fav_VH_feature_names",row.names=F,col.names=F,sep = ',')

##### variable selection with "MD" or minimal depth method #####
DFS_Clin_fav_MD_features<- var.select.rfsrc(Surv(DFI.time, DFI) ~ ., DFS_Clin_fav,method = "md", nrep = 100, nstep = 5)

#### Extract names of variables #####

DFS_Clin_fav_MD_feature_names = DFS_Clin_fav_MD_features$topvars

write.table(DFS_Clin_fav_MD_feature_names,file="DFS_Clin_fav_MD_feature_names",row.names=F,col.names=F,sep = ',')



##### variable selection with "vh.vimp" or variable hunting with VIMP (variable importance). method #####
DFS_Clin_fav_vh.vimp_features<- var.select.rfsrc(Surv(DFI.time, DFI) ~ ., DFS_Clin_fav,method = "vh.vimp", nrep = 100, nstep = 5)

#### Extract names of variables #####

DFS_Clin_fav_VH_vimp_feature_names = DFS_Clin_fav_vh.vimp_features$topvars

write.table(DFS_Clin_fav_VH_vimp_feature_names,file="DFS_Clin_fav_VH_vimp_feature_names",row.names=F,col.names=F,sep = ',')


###########################################################################################################################

######################## feature selection for unfavorable features ####################################

##### variable selection with "VH" or variable hunting method #####
DFS_Clin_unfav_VH_features<- var.select.rfsrc(Surv(DFI.time, DFI) ~ ., DFS_Clin_unfav,method = "vh", nrep = 100, nstep = 5)

#### Extract names of variables #####

DFS_Clin_unfav_VH_feature_names <- DFS_Clin_unfav_VH_features$topvars

write.table(DFS_Clin_unfav_VH_feature_names,file="DFS_Clin_unfav_VH_feature_names",row.names=F,col.names=F,sep = ',')

##### variable selection with "MD" or minimal depth method #####
DFS_Clin_unfav_MD_features<- var.select.rfsrc(Surv(DFI.time, DFI) ~ ., DFS_Clin_unfav,method = "md", nrep = 100, nstep = 5)

#### Extract names of variables #####

DFS_Clin_unfav_MD_feature_names = DFS_Clin_unfav_MD_features$topvars

write.table(DFS_Clin_unfav_MD_feature_names,file="DFS_Clin_unfav_MD_feature_names",row.names=F,col.names=F,sep = ',')



##### variable selection with "vh.vimp" or variable hunting with VIMP (variable importance). method #####
DFS_Clin_unfav_vh.vimp_features<- var.select.rfsrc(Surv(DFI.time, DFI) ~ ., DFS_Clin_unfav,method = "vh.vimp", nrep = 100, nstep = 5)

#### Extract names of variables #####

DFS_Clin_unfav_VH_vimp_feature_names = DFS_Clin_unfav_vh.vimp_features$topvars

write.table(DFS_Clin_unfav_VH_vimp_feature_names,file="DFS_Clin_unfav_VH_vimp_feature_names",row.names=F,col.names=F,sep = ',')


#############################################################################################################################



######################## feature selection for combined fav and unfav features ####################################

##### variable selection with "VH" or variable hunting method #####
DFS_Clin_both_fav_unfav_VH_features<- var.select.rfsrc(Surv(DFI.time, DFI) ~ ., DFS_Clin_both_fav_unfav,method = "vh", nrep = 100, nstep = 5)

#### Extract names of variables #####

DFS_Clin_both_VH_feature_names <- DFS_Clin_both_fav_unfav_VH_features$topvars

write.table(DFS_Clin_both_VH_feature_names,file="DFS_Clin_both_VH_feature_names",row.names=F,col.names=F,sep = ',')

##### variable selection with "MD" or minimal depth method #####
DFS_Clin_both_fav_unfav_MD_features<- var.select.rfsrc(Surv(DFI.time, DFI) ~ ., DFS_Clin_both_fav_unfav,method = "md", nrep = 100, nstep = 5)

#### Extract names of variables #####

DFS_Clin_both_MD_feature_names = DFS_Clin_both_fav_unfav_MD_features$topvars

write.table(DFS_Clin_both_MD_feature_names,file="DFS_Clin_both_MD_feature_names",row.names=F,col.names=F,sep = ',')



##### variable selection with "vh.vimp" or variable hunting with VIMP (variable importance). method #####
DFS_Clin_both_fav_unfav_vh.vimp_features<- var.select.rfsrc(Surv(DFI.time, DFI) ~ ., DFS_Clin_both_fav_unfav,method = "vh.vimp", nrep = 100, nstep = 5)

#### Extract names of variables #####

DFS_Clin_both_VH_vimp_feature_names = DFS_Clin_both_fav_unfav_vh.vimp_features$topvars

write.table(DFS_Clin_both_VH_vimp_feature_names,file="DFS_Clin_both_VH_vimp_feature_names",row.names=F,col.names=F,sep = ',')


