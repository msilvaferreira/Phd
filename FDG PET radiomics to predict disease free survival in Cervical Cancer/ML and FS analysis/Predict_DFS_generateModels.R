rm(list = ls())


time1 <- proc.time()

## Libraries##

library(prodlim)
library(Publish)
library(caret)
library(glmnet)
library(pROC)
library(DMwR)
library(randomForest)
library(modelr)
library(MLmetrics)
library(boot)
library(praznik)
library(svDialogs)
library(e1071)
library(data.table)
library(ROCR)
library(ROSE)
library(matrixStats)


General_path='C:\\Users\\f\\Documents\\Phd\\Project 1\\'

########## Functions###################
source(paste(General_path,'code\\Paper 1\\Load_data.R',sep=""))
source(paste(General_path,'code\\Paper 1\\Transform_ComBat.R',sep=""))
source(paste(General_path,'code\\Paper 1\\CV_Select_FS_NrFeatures.R',sep=""))
source(paste(General_path,'code\\Paper 1\\Test_5CV_FS_NrFeatures.R',sep=""))
source(paste(General_path,'code\\Paper 1\\MRMR.R',sep=""))
source(paste(General_path,'code\\Paper 1\\FeatureSelection_RandomForest.R',sep=""))
source(paste(General_path,'code\\Paper 1\\my_RF_to_BOOTSTRAP_Sensitivity_AUC.R',sep=""))
source(paste(General_path,'code\\Paper 1\\Discretise.R',sep=""))


#######################################################################

####### Define Task with propmpt/ answer system ##########

#cat("Which data you want to use as test data?\n1)Mix data sets from different institutes\n2)A\n3)B\n4)EF\n5)G\6)H")
Var2=c(6)
#cat("Which classifier you want to use?\n1)RF\n2)SVM\n3)Logistic Regression\n4)Naive Bayes")
Var3=c(1)

opt_criteria='Mean.fscore'
beta=1
kernel_svm='radial'


#cat("Which feature selection method do you want to use?\n1)Forward MRMR pearson\n2)Backward MRMR pearson\n3)Forward MRMR spearman\n4)Backward MRMR spearman\n5)Forward MRMR MI\n6)Embebed RF FS Accuracy decrease\n7)Embebed RF FS Gini Impurity\n8)kmeans_corr\n9)PCA\n")
Var4=c(5)

#cat("Which features do you want to use?\n1)Tumor Radiomics\n2)Tumor to liver Ratio radiomics\n")
Var5=c(2)
#cat("Do you want to add any other information?\n1)Only Radiomics\n2)Radiomics + clinical\n3)Radiomics with interpolated images\n4) Radiomics with Combat correction\n")
Var6=c(2)

#cat("Which discretisation method you want to use?\n1)Fix bin width\n2)Fix bin number\n")
Var7=c(1)
if ((1 %in% Var7)==TRUE){
  #cat("Which discretisation width?\n1)0.5\n2)0.2\n3)0.1\n4)0.05\n5)All\n")
  Var8=c(4)
}
if ((2 %in% Var7)==TRUE){
  #cat("Which discretisation number?\n1)64\n2)32\n3)All\n")
  Var9=c(2)
}




c_df=1



for (i2 in 1:length(Var2)) {
  for (i3 in 1:length(Var3)) {
    for (i4 in 1:length(Var4)) {
      for (i5 in 1:length(Var5)) {
        for (i6 in 1:length(Var6)) {
          for (i7 in 1:length(Var7)) {
            
            if (Var7[i7]==1 ){
              Var89=Var8
            }
            if (Var7[i7]==2){
              Var89=Var9
            }
            for (i89 in 1:length(Var89)) {
              
              
              #######################################################################
              
              ##############################################
              ############Load Data ########################
              Data_Patient= read.csv(paste(General_path,'Excell_data\\Other\\Data_169patient.csv',sep=""), sep = ';')
              Data_Patient= data.frame("Patient" = Data_Patient$Patient)
              
              DATA_C=Load_data(General_path,"clinical", Data_Patient, 'None','None')
              DATA_C=merge(DATA_C,Data_Patient, by='Patient')
              DATA_C$Institute=NULL
              
              DATA_output=Load_data(General_path,"Output", Data_Patient, 'None','None')
              DATA_output=merge(DATA_output,Data_Patient, by='Patient')
              DATA_output_DFS=data.frame('Patient'=DATA_output$Patient, 'Outcome'=DATA_output$DFS)
              
              Data_scanner=Load_data(General_path,"Scanner", Data_Patient, 'None','None')
              Data_scanner=merge(Data_scanner,Data_Patient,by='Patient')
              
              
              if (Var7[i7]==1){
                Disccretisation_method='FBW'
              }
              if (Var7[i7]==2){
                Disccretisation_method='FBN'
              }
              
              if (Var6[i6]==3){
                interpolation=TRUE
              }else{
                interpolation=FALSE
              }
              
              if (interpolation==FALSE){
                if (Disccretisation_method=='FBW'){
                  if (Var5[i5]==2){
                    DATA_radiomics_all=Load_data(General_path,"TLR_Radiomics",Data_Patient, 'FBW','False')
                    
                  }else{
                    DATA_radiomics_all=Load_data(General_path,"OR_Radiomics", Data_Patient, 'FBW','False')
                  }
                }
                if (Disccretisation_method=='FBN'){
                  if (Var5[i5]==2){
                    DATA_radiomics_all=Load_data(General_path,"TLR_Radiomics",Data_Patient, 'FBN','False')
                  } else{
                    DATA_radiomics_all=Load_data(General_path,"OR_Radiomics", Data_Patient, 'FBN','False')
                  }
                }
              }
              if (interpolation==TRUE){
                if (Disccretisation_method=='FBW'){
                  if (Var5[i5]==2){
                    DATA_radiomics_all=Load_data(General_path,"TLR_Radiomics",Data_Patient, 'FBW','True')
                    
                  }else{
                    DATA_radiomics_all=Load_data(General_path,"OR_Radiomics", Data_Patient, 'FBW','True')
                  }
                }
                if (Disccretisation_method=='FBN'){
                  if (Var5[i5]==2){
                    DATA_radiomics_all=Load_data(General_path,"TLR_Radiomics",Data_Patient, 'FBN','True')
                  } else{
                    DATA_radiomics_all=Load_data(General_path,"OR_Radiomics", Data_Patient, 'FBN','True')
                  }
                }
              }
              
              
              
              
              #####################################################################################
              ######################### Choose discretisation #####################################
              
              if (Var7[i7]==1){
                if (Var89[i89]==1){
                  Disccretisation=0.5
                }
                if (Var89[i89]==2){
                  Disccretisation=0.2
                }
                if (Var89[i89]==3){
                  Disccretisation=0.1
                }
                if (Var89[i89]==4){
                  Disccretisation=0.05
                }
                if (Var89[i89]==5){
                  Disccretisation='All'
                }
              }
              
              if (Var7[i7]==2){
                if (Var89[i89]==1){
                  Disccretisation=64
                }
                if (Var89[i89]==2){
                  Disccretisation=32
                }
                if (Var89[i89]==3){
                  Disccretisation='All'
                }
              }
              
              DATA_radiomics_all=Discretise(DATA_radiomics_all,Disccretisation)
              
              dp=DATA_radiomics_all$Patient
              DATA_radiomics_all$Patient=NULL
              DATA_radiomics_all=DATA_radiomics_all[is.finite(colSums(DATA_radiomics_all)) ]
              
              DATA_radiomics_all =DATA_radiomics_all[, sapply(data.frame(sapply(DATA_radiomics_all, as.numeric)), var) != 0]
              
              
              
              
              #####################################################################################
              
              ###########combat harmonization#################
              
              if (Var6[i6]==4 ){
                
                DATA_radiomics_all=Transform_ComBat(DATA_radiomics_all, Data_Patient, Data_scanner)
              }
              
              
              
              #####################################################################################
              
              #########################################################################
              ########## Merge radiomics with clinical Output###########################
              
              
              DATA_radiomics_all_Patient=DATA_radiomics_all
              DATA_radiomics_all_Patient['Patient']=dp
              df_inter=merge(Data_Patient,DATA_output,by="Patient")
              DATA_radiomics_plus_Outcome<- merge( data.frame('Patient'=df_inter$Patient,'Outcome'=df_inter$DFS),DATA_radiomics_all_Patient,by="Patient")
              
              
              
              
              
              ##################################################################################
              
              ##################################################################################
              ###########################Adding clinical data ##################################
              
              
              if (Var6[i6]==2){  
                DATA_C$Institute <- NULL
                DATA_radiomics_plus_Outcome<- merge( DATA_C,DATA_radiomics_plus_Outcome,by="Patient")
              }
              
              
              ###############################################################################
              
              #####################################################################
              #######################Data to be used###############################
              
              DATA_radiomics_plus_Outcome$X <- NULL
              DATA_radiomics_plus_Outcome['Outcome']=as.factor(DATA_radiomics_plus_Outcome$Outcome)
              Data2= DATA_radiomics_plus_Outcome
              L_P_Data=Data2$Patient
              DATA_radiomics_plus_Outcome$Patient <- NULL
              DATA_radiomics_plus_Outcome=DATA_radiomics_plus_Outcome[is.finite(colSums(DATA_radiomics_plus_Outcome[, names(DATA_radiomics_plus_Outcome) != "Outcome"])) ]
              DATA_radiomics_plus_Outcome=DATA_radiomics_plus_Outcome[, sapply(data.frame(sapply(DATA_radiomics_plus_Outcome, as.numeric)), var) != 0]
              DATA_radiomics_plus_Outcome['Patient']=L_P_Data
              Data= DATA_radiomics_plus_Outcome
              
              
              
              #Data$Patient <- NULL
              Data$Month <- NULL
              
              Data_R_F= Data
              Z=1
              Rep=1
              AUC_general_cross <- matrix(0, Rep, 1)
              AUC_general_test <- matrix(0, Rep, 1)
              AUC_general_train <- matrix(0, Rep, 1)
              cc=1
              
              
              
              
              
              ###############################################################################
              
              ###############################################################################
              #######################Split train/test/val####################################
              
              
              
              D_A=Data_scanner[which(Data_scanner$Scanner=='A'),]
              
              D_EF=Data_scanner[which(Data_scanner$Scanner=='EF'),]
              
              D_H=Data_scanner[which(Data_scanner$Scanner=='H'),]
              
              ratio=0.8
              
              set.seed(Z)
              train_A <- sample(nrow(D_A), ratio*nrow(D_A), replace = FALSE)#
              Set_A <- D_A[train_A,]
              TestSet_A <- D_A[-train_A,]
              
              
              set.seed(Z)
              train_EF <- sample(nrow(D_EF), ratio*nrow(D_EF), replace = FALSE)#
              Set_EF <- D_EF[train_EF,]
              TestSet_EF <- D_EF[-train_EF,]
              
              TestSet_H=D_H
              
              Set=rbind(Set_A,Set_EF)
              Set$Institute=NULL
              Set$Scanner=NULL
              
              Set=merge(Set,Data_R_F,by='Patient')
              
              Set$Patient <- NULL
              
              Set=rbind(Set,Set[which(Set$Outcome==1),])
              
              if(Var2[i2]==1){
                
                TestSet=rbind(TestSet_A,TestSet_EF)
                TestSet$Institute=NULL
                TestSet$Scanner=NULL
                
                TestSet=merge(TestSet,Data_R_F,by='Patient')
                TestSet$Patient <- NULL
                
                
                
              }else{
                
                if (Var2[i2]==2){
                  TestSet=TestSet_A
                  TestSet$Institute=NULL
                  TestSet$Scanner=NULL
                }
                
                if (Var2[i2]==4){
                  TestSet=TestSet_EF
                  TestSet$Institute=NULL
                  TestSet$Scanner=NULL
                }
                
                if (Var2[i2]==6){
                  TestSet=TestSet_H
                  TestSet$Institute=NULL
                  TestSet$Scanner=NULL
                }
                
                
                TestSet=merge(TestSet,Data_R_F,by='Patient')
                TestSet$Patient <- NULL
                
                
                
              }
              
              
              #############################################################################################
              
              ##############################################################################################
              #######################5 CV for tunning number of features####################################
              
              
              folds <- 5
              v_f=seq(1, 10 )
              N <- length(v_f)
              
              if (Var4[i4]==1){
                direction='up'
                FS='MRMR_pearson'#MRMR_pearson, RF, MRMR_MI, MRMR_spearman
                criteria='0'
              }
              if (Var4[i4]==2){
                direction='down'
                FS='MRMR_pearson'#MRMR_pearson, RF, MRMR_MI, MRMR_spearman
                criteria='0'
              }
              if (Var4[i4]==3){
                direction='up'
                FS='MRMR_spearman'#MRMR_pearson, RF, MRMR_MI, MRMR_spearman
                criteria='0'
              }
              if (Var4[i4]==4){
                direction='down'
                FS='MRMR_spearman'#MRMR_pearson, RF, MRMR_MI, MRMR_spearman
                criteria='0'
              }
              if (Var4[i4]==5){
                direction='up'
                FS='MRMR_MI'#MRMR_pearson, RF, MRMR_MI, MRMR_spearman
                criteria='0'
              }
              if (Var4[i4]==6){
                criteria='MeanDecreaseAccuracy'#'%IncMSE' 'IncNodePurity'
                FS='RF'#MRMR_pearson, RF, MRMR_MI, MRMR_spearman
              }
              if (Var4[i4]==7){
                criteria='MeanDecreaseGini'#'%IncMSE' 'IncNodePurity'
                FS='RF'#MRMR_pearson, RF, MRMR_MI, MRMR_spearman
              }
              
              
              
              set.seed(Z)
              folding <- crossv_kfold(Set,k=folds)
              
              
              Res=CV_Select_FS_NrFeatures(Z,Set,N,folds, FS, Var3[i3],folding, criteria,opt_criteria,v_f,beta)
              
              NR_Features=Res[[1]]
              
              C1_DF_10_auc_order=Res[opt_criteria][[1]] #Mean.AUC
              
              
              
              
              
              ##############################################################################################
              #######################Generate desired features############################################
              
              
              
              TR_Features=Set
              TR_Features_out=TR_Features['Outcome']
              TR_Features['Outcome']=NULL
              TR_Features =TR_Features[, sapply(data.frame(sapply(TR_Features, as.numeric)), var) != 0]
              
              
              if ((FS=='MRMR_pearson')|(FS=='MRMR_spearman')){
                TR_Features=FeatureSelection_MLMR(TR_Features,Set$Outcome,NR_Features,direction,FS)
              }
              if (FS=='RF'){
                TR_Features=FeatureSelection_RandomForest(TR_Features,Z,Set$Outcome,criteria,NR_Features)
              }
              
              if (FS=='MRMR_MI'){
                lid <- matrix(0, NR_Features, 1)
                old_f=colnames(TR_Features)
                TR_F=MRMR(TR_Features, Set$Outcome, k = NR_Features)
                for (nn in 1:(length(TR_F$selection))){
                  lid[nn]=TR_F$selection[[nn]]
                }
                features=old_f[lid]
                TR_Features=TR_Features[c(features)]
              }
              
              
              features=colnames(TR_Features)
              print(features)
              
              
              
              
              
              Res_Test=Test_5CV_FS_NrFeatures(TestSet, features,Set, folding, folds, Var3[i3], Z)
              CI_down_TEST=Res_Test[[1]]
              CI_up_TEST= Res_Test[[2]]
              mean_TEST=Res_Test[[3]]
              
              CI_down_TEST_fscore=Res_Test[[4]]
              CI_up_TEST_fscore= Res_Test[[5]]
              mean_TEST_fscore=Res_Test[[6]]
              mean_TEST_precision=Res_Test[[7]]
              mean_TEST_recall=Res_Test[[8]]
              
              ##################################################
              ########## Test with bootstrap############
              
              
              
              TrainSet_out= Set$Outcome
              
              SetSet=Set
              
              SetSet=SetSet[, sapply(data.frame(sapply(SetSet, as.numeric)), var) != 0]
              
              TrainSet_new= SetSet[c(features)]
              
              
              TrainSet_new['Outcome']=TrainSet_out
              
              TestSet_out= TestSet$Outcome
              
              TestSetTestSet=TestSet
              
              
              TestSet_new= TestSetTestSet[c(features)]
              
              TestSet_new['Outcome']=TestSet_out
              
              
              
              
              
              my_bootstrap=boot(TrainSet_new,my_RF_to_BOOTSTRAP_Sensitivity_AUC, R=1000, TEST=TestSet_new,Var3=Var3[i3], features=features,beta=beta)
              
              auc_train_mean=mean(my_bootstrap$t[,1])
              auc_test_mean=mean(my_bootstrap$t[,2])
              fscore_train_mean=mean(my_bootstrap$t[,10])
              fscore_test_mean=mean(my_bootstrap$t[,11])
              prauc_test_mean=mean(my_bootstrap$t[,14])
              precision_test_mean=mean(my_bootstrap$t[,12])
              recall_test_mean=mean(my_bootstrap$t[,13])
              #sensitivity_train_mean=mean(my_bootstrap$t[,3])
              #sensitivity_test_mean=mean(my_bootstrap$t[,4])
              #specificity_train_mean=mean(my_bootstrap$t[,5])
              #specificity_test_mean=mean(my_bootstrap$t[,6])
              #Accuracy_test_mean=mean(my_bootstrap$t[,8])
              
              #CI_TRAIN=boot.ci(my_bootstrap, index=1)
              
              if(class(try(boot.ci(my_bootstrap, index=2), silent=FALSE))=="try-error"|| class(try(boot.ci(my_bootstrap, index=11), silent=FALSE))=="try-error")
              {
                
              }else{
                CI_TEST=boot.ci(my_bootstrap, index=2)
                CI_TEST_fscore=boot.ci(my_bootstrap, index=11)
                CI_TEST_prauc=boot.ci(my_bootstrap, index=14)
              }
              #CI_sensitivity_TRAIN=boot.ci(my_bootstrap, index=3)
              #CI_sensitivity_TEST=boot.ci(my_bootstrap, index=4)
              #CI_specificity_TRAIN=boot.ci(my_bootstrap, index=5)
              #CI_specificity_TEST=boot.ci(my_bootstrap, index=6)
              #CI_Accuracy_TEST=boot.ci(my_bootstrap, index=8)
              
              
              
              
              ############################################################################################
              
              ############################################################################################
              #######################Save results in data frame############################################
              
              if (Var2[i2]==1){
                Var2_df='Mix'
              }
              if (Var2[i2]==2){
                Var2_df='A'
              }
              if (Var2[i2]==3){
                Var2_df='B'
              }
              if (Var2[i2]==4){
                Var2_df='EF'
              }
              if (Var2[i2]==5){
                Var2_df='G'
              }
              if (Var2[i2]==6){
                Var2_df='H'
              }
              
              if (Var3[i3]==1){
                Var3_df='RF'
              }
              if (Var3[i3]==2){
                Var3_df='SVM'
              }
              if (Var3[i3]==3){
                Var3_df='Logistic Regression'
              }
              if (Var3[i3]==4){
                Var3_df='Naive Bayes'
              }
              
              if (Var4[i4]==1){
                Var4_df='Forward MRMR pearson'
              }
              if (Var4[i4]==2){
                Var4_df='Backward MRMR pearson'
              }
              if (Var4[i4]==3){
                Var4_df='Forward MRMR spearman'
              }
              if (Var4[i4]==4){
                Var4_df='Backward MRMR spearman'
              }
              if (Var4[i4]==5){
                Var4_df='Forward MRMR MI'
              }
              if (Var4[i4]==6){
                Var4_df='Embebed RF FS Accuracy decrease'
              }
              if (Var4[i4]==7){
                Var4_df='Embebed RF FS Gini Impurity'
              }
              if (Var4[i4]==8){
                Var4_df='kmeans_corr'
              }
              if (Var4[i4]==9){
                Var4_df='PCA'
              }
              
              if (Var5[i5]==1){
                Var5_df='Original Radiomics Features'
              }
              if (Var5[i5]==2){
                Var5_df='Tumor to Liver ratio Radiomics Features'
              }
              
              if (Var6[i6]==1){
                Var6_df='Empty'
              }
              if (Var6[i6]==2){
                Var6_df='Clinical'
              }
              if (Var6[i6]==3){
                Var6_df='Interpolation'
              }
              if (Var6[i6]==4){
                Var6_df='Combat harmonization'
              }
              
              if (Var7[i7]==1){
                Var7_df='FBW'
                if (Var89[i89]==1){
                  Var89_df=0.5
                }
                if (Var89[i89]==2){
                  Var89_df=0.2
                }
                if (Var89[i89]==3){
                  Var89_df=0.1
                }
                if (Var89[i89]==4){
                  Var89_df=0.05
                }
                if (Var89[i89]==5){
                  Var89_df='All'
                }
              }
              if (Var7[i7]==2){
                Var7_df='FBN'
                if (Var89[i89]==1){
                  Var89_df=64
                }
                if (Var89[i89]==2){
                  Var89_df=32
                }
                if (Var89[i89]==3){
                  Var89_df='All'
                }
              }
              
              CV_Val_fscore=paste(toString(round(Res$Mean.fscore[length(C1_DF_10_auc_order)],2)),'(',toString(round(Res$CI_down_fscore[length(C1_DF_10_auc_order)],2)),'-',toString(round(Res$CI_up_fscore[length(C1_DF_10_auc_order)],2)),')', sep = " ", collapse = NULL)
              CVR5_fscore=paste(toString(round(mean_TEST_fscore,2)),'(',toString(round(CI_down_TEST_fscore,2)),'-',toString(round(CI_up_TEST_fscore,2)),')', sep = " ", collapse = NULL)
              if (mean_TEST_fscore>=0){
                
                error_p=try(paste(toString(round(prauc_test_mean,2),round(fscore_test_mean,2)),'(',toString(round((CI_TEST_fscore[6])$percent[4],2)),'-',toString(round((CI_TEST_fscore[6])$percent[5],2)),')', sep = " ", collapse = NULL),silent=T)
                if (class(error_p)=="try-error"){
                  Boot_test_fscore=paste(toString(round(fscore_test_mean,2)),'(',toString('Empty'),'-',toString('Empty'),')', sep = " ", collapse = NULL)
                  Boot_test_prauc=paste(toString(round(prauc_test_mean,2)),'(',toString('Empty'),'-',toString('Empty'),')', sep = " ", collapse = NULL)
                  
                  B_min_fscore='Empty'
                  B_max_fscore='Empty'
                  B_mean_fscore=round(fscore_test_mean,2)
                  B_mean_prauc=round(prauc_test_mean,2)
                  B_mean_prauc=round(prauc_test_mean,2)
                  B_mean_test_precision=round(precision_test_mean,2)
                  B_mean_test_recall=round(recall_test_mean,2)
                  
                  
                  
                }else{
                  Boot_test_fscore=paste(toString(round(fscore_test_mean,2)),'(',toString(round((CI_TEST_fscore[6])$percent[4],2)),'-',toString(round((CI_TEST_fscore[6])$percent[5],2)),')', sep = " ", collapse = NULL)
                  Boot_test_prauc=toString(round(prauc_test_mean,2))
                  
                  B_min_fscore=round((CI_TEST_fscore[6])$percent[4],2)
                  B_max_fscore=round((CI_TEST_fscore[6])$percent[5],2)
                  B_mean_fscore=round(fscore_test_mean,2)
                  B_mean_prauc=round(prauc_test_mean,2)
                  B_mean_test_precision=round(precision_test_mean,2)
                  B_mean_test_recall=round(recall_test_mean,2)
                  
                  
                  
                }
              }else{
                Boot_test_fscore='Empty'
                B_min_fscore='Empty'
                B_max_fscore='Empty'
                B_mean_fscore='Empty'
                B_mean_prauc='Empty'
                B_mean_test_precision='Empty'
                B_mean_test_recall='Empty'
                
                
              }
              
              CV_Val=paste(toString(round(Res$Mean.AUC[length(C1_DF_10_auc_order)],2)),'(',toString(round(Res$CI_down[length(C1_DF_10_auc_order)],2)),'-',toString(round(Res$CI_up[length(C1_DF_10_auc_order)],2)),')', sep = " ", collapse = NULL)
              CVR5=paste(toString(round(mean_TEST,2)),'(',toString(round(CI_down_TEST,2)),'-',toString(round(CI_up_TEST,2)),')', sep = " ", collapse = NULL)
              
              
              if (mean_TEST>=0){
                
                error_p=try(paste(toString(round(auc_test_mean,2)),'(',toString(round((CI_TEST[6])$percent[4],2)),'-',toString(round((CI_TEST[6])$percent[5],2)),')', sep = " ", collapse = NULL),silent=T)
                if (class(error_p)=="try-error"){
                  Boot_test=paste(toString(round(auc_test_mean,2)),'(',toString('Empty'),'-',toString('Empty'),')', sep = " ", collapse = NULL)
                  B_min='Empty'
                  B_max='Empty'
                  B_mean=round(auc_test_mean,2)
                  
                  
                  
                }else{
                  Boot_test=paste(toString(round(auc_test_mean,2)),'(',toString(round((CI_TEST[6])$percent[4],2)),'-',toString(round((CI_TEST[6])$percent[5],2)),')', sep = " ", collapse = NULL)
                  B_min=round((CI_TEST[6])$percent[4],2)
                  B_max=round((CI_TEST[6])$percent[5],2)
                  B_mean=round(auc_test_mean,2)
                  
                  
                  
                }
              }else{
                Boot_test='Empty'
                B_min='Empty'
                B_max='Empty'
                B_mean='Empty'
                
                
                
              }
              
              print(paste(toString(Var2_df),'_',toString(Var3_df),'_',toString(Var4_df),'_',toString(Var5_df),'_',toString(Var6_df),'_',toString(Var7_df),'_',toString(Var89_df), sep = " ", collapse = NULL))
              B_mean_train=auc_train_mean
              B_mean_train_fscore=fscore_train_mean
              if (c_df==1){
                Results_df=data.frame('Test data'=Var2_df, 'Classifier'= Var3_df, 'Feature Selection'=Var4_df, 'Features'=Var5_df, 'Additional Information'=Var6_df, 'Discretisation method'=Var7_df, 'Discretisation width/number'=Var89_df,'Number of features'=NR_Features,'CV5_Val_AUC'=CV_Val,'CV5_Test_AUC'=CVR5,'Boot_Test_AUC'=Boot_test,'Boot train'=B_mean_train, 'Features'= toString(features))
                Results_df2=data.frame('Test data'=Var2_df, 'Classifier'= Var3_df, 'Feature Selection'=Var4_df, 'Features'=Var5_df, 'Additional Information'=Var6_df, 'Discretisation method'=Var7_df, 'Discretisation width/number'=Var89_df,'Number of features'=NR_Features,'CV5_Val_AUC_mean'=round(Res$Mean.AUC[length(C1_DF_10_auc_order)],2),'CV5_Val_AUC_min'=round(Res$CI_down[length(C1_DF_10_auc_order)],2),'CV5_Val_AUC_max'=round(Res$CI_up[length(C1_DF_10_auc_order)],2),'CV5_Test_AUC'=round(mean_TEST,2),'CV5_Test_AUC_min'=round(CI_down_TEST,2),'CV5_Test_AUC_max'=round(CI_up_TEST,2),'Boot_Test_AUC'=B_mean, 'Boot_Test_AUC_min'=B_min,'Boot_Test_AUC_max'=B_max,'Boot train'=B_mean_train, 'Boot_Test_AUC_fscore'=toString(B_mean_fscore),'Boot_Test_AUC_min_fscore'=B_min_fscore,'Boot_Test_AUC_max_fscore'=B_max_fscore,'Boot train fscore'=B_mean_train_fscore, 'Boot_Test_precision_fscore'=B_mean_test_precision, 'Boot_Test_recall_fscore'=B_mean_test_recall ,'Boot_Test_AUC_prauc'=toString(B_mean_prauc), 'Features'= toString(features))
                Results_df3=data.frame('Test data'=Var2_df, 'Classifier'= Var3_df, 'Feature Selection'=Var4_df, 'Features'=Var5_df, 'Additional Information'=Var6_df, 'Discretisation method'=Var7_df, 'Discretisation width/number'=Var89_df,'Number of features'=NR_Features,'CV5_Val_fscore_mean'=round(Res$Mean.fscore[length(C1_DF_10_auc_order)],2),'CV5_Val_fscore_min'=round(Res$CI_down_fscore[length(C1_DF_10_auc_order)],2),'CV5_Val_fscore_max'=round(Res$CI_up_fscore[length(C1_DF_10_auc_order)],2),'CV5_Test_fscore'=round(mean_TEST_fscore,2),'CV5_Test_fscore_min'=round(CI_down_TEST_fscore,2),'CV5_Test_fscore_max'=round(CI_up_TEST_fscore,2), 'CV5_Test_precision'=toString(round(mean_TEST_precision,2)),'CV5_Test_recall'=toString(round(mean_TEST_recall,2)),'Boot_Test_AUC_fscore'=toString(B_mean_fscore),'Boot_Test_AUC_min_fscore'=B_min_fscore,'Boot_Test_AUC_max_fscore'=B_max_fscore,'Boot train fscore'=B_mean_train_fscore, 'Boot_Test_precision_fscore'=B_mean_test_precision ,'Boot_Test_recall_fscore'=B_mean_test_recall , 'Boot_Test_AUC_prauc'=toString(B_mean_prauc), 'Features'= toString(features))
                
                print(B_mean)
                c_df=c_df+1
              }else{
                Results_df_new=data.frame('Test data'=Var2_df, 'Classifier'= Var3_df, 'Feature Selection'=Var4_df, 'Features'=Var5_df, 'Additional Information'=Var6_df, 'Discretisation method'=Var7_df, 'Discretisation width/number'=Var89_df,'Number of features'=NR_Features,'CV5_Val_AUC'=CV_Val,'CV5_Test_AUC'=CVR5,'Boot_Test_AUC'=Boot_test,'Boot train'=B_mean_train, 'Features'= toString(features))
                Results_df <- rbind(Results_df, Results_df_new)
                
                Results_df2_new=data.frame('Test data'=Var2_df, 'Classifier'= Var3_df, 'Feature Selection'=Var4_df, 'Features'=Var5_df, 'Additional Information'=Var6_df, 'Discretisation method'=Var7_df, 'Discretisation width/number'=Var89_df,'Number of features'=NR_Features,'CV5_Val_AUC_mean'=round(Res$Mean.AUC[length(C1_DF_10_auc_order)],2),'CV5_Val_AUC_min'=round(Res$CI_down[length(C1_DF_10_auc_order)],2),'CV5_Val_AUC_max'=round(Res$CI_up[length(C1_DF_10_auc_order)],2),'CV5_Test_AUC'=toString(round(mean_TEST,2)),'CV5_Test_AUC_min'=toString(round(CI_down_TEST,2)),'CV5_Test_AUC_max'=toString(round(CI_up_TEST,2)),'Boot_Test_AUC'=toString(B_mean), 'Boot_Test_AUC_min'=toString(B_min),'Boot_Test_AUC_max'=toString(B_max),'Boot train'=B_mean_train, 'Boot_Test_AUC_fscore'=toString(B_mean_fscore),'Boot_Test_AUC_min_fscore'=toString(B_min_fscore),'Boot_Test_AUC_max_fscore'=toString(B_max_fscore),'Boot train fscore'=B_mean_train_fscore,'Boot_Test_precision_fscore'=B_mean_test_precision ,'Boot_Test_recall_fscore'=B_mean_test_recall ,'Boot_Test_AUC_prauc'=toString(B_mean_prauc),'Features'= toString(features))
                Results_df2 <- rbind(Results_df2, Results_df2_new)
                
                Results_df3_new=data.frame('Test data'=Var2_df, 'Classifier'= Var3_df, 'Feature Selection'=Var4_df, 'Features'=Var5_df, 'Additional Information'=Var6_df, 'Discretisation method'=Var7_df, 'Discretisation width/number'=Var89_df,'Number of features'=NR_Features,'CV5_Val_fscore_mean'=round(Res$Mean.fscore[length(C1_DF_10_auc_order)],2),'CV5_Val_fscore_min'=round(Res$CI_down_fscore[length(C1_DF_10_auc_order)],2),'CV5_Val_fscore_max'=round(Res$CI_up_fscore[length(C1_DF_10_auc_order)],2),'CV5_Test_fscore'=toString(round(mean_TEST_fscore,2)),'CV5_Test_fscore_min'=toString(round(CI_down_TEST_fscore,2)),'CV5_Test_fscore_max'=toString(round(CI_up_TEST_fscore,2)),'CV5_Test_precision'=toString(round(mean_TEST_precision,2)),'CV5_Test_recall'=toString(round(mean_TEST_recall,2)), 'Boot_Test_AUC_fscore'=toString(B_mean_fscore),'Boot_Test_AUC_min_fscore'=toString(B_min_fscore),'Boot_Test_AUC_max_fscore'=toString(B_max_fscore),'Boot train fscore'=B_mean_train_fscore,'Boot_Test_precision_fscore'=B_mean_test_precision ,'Boot_Test_recall_fscore'=B_mean_test_recall, 'Boot_Test_AUC_prauc'=toString(B_mean_prauc), 'Features'= toString(features))
                Results_df3 <- rbind(Results_df3, Results_df3_new)
                
                print(B_mean)
                print(B_mean_fscore)
              }
            }
            
            
          }
        }
      }
    }
  }
}


time= proc.time() - time1

# write.csv(Results_df2, file = "Test1.csv")


