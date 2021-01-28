rm(list = ls())

##### libraries ######

library(survival)
library(survminer)
library(standardize)
library(e1071)
library(ROCR)
library(boot)
library(pRROC)
library(caret)

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

### FUNCTIONS###


##### load functions #####
General_path='C:\\Users\\f\\Documents\\Phd\\Project 1\\'


source(paste(General_path,'code\\Paper2\\Load_data.R',sep=""))
source(paste(General_path,'code\\Paper2\\Performance_TestData.R',sep=""))
source(paste(General_path,'code\\Paper2\\Performance_class.R',sep=""))
source(paste(General_path,'code\\Paper2\\Performance_Total.R',sep=""))
source(paste(General_path,'code\\Paper2\\chose_descretisation_nr.R',sep=""))
source(paste(General_path,'code\\Paper 1\\MRMR.R',sep=""))
source(paste(General_path,'code\\Paper 1\\FeatureSelection_RandomForest.R',sep=""))




#### read data ######


Data_Patient= read.csv("C:\\Users\\f\\Documents\\Phd\\Project 1\\Excell_data\\Other\\Data_169patient.csv", sep = ';')
Data_Patient= data.frame("Patient" = Data_Patient$Patient)

Data_scanner=read.csv("C:\\Users\\f\\Documents\\Phd\\Project 1\\Excell_data\\Other\\Data169_scanner_binary.csv", sep = ',')
Data_scanner=merge(Data_scanner,Data_Patient,by='Patient')

DATA_C=Load_data(General_path,"clinical", Data_Patient, 'None','None')
DATA_C$Institute=NULL

DATA_output=Load_data(General_path,"Output", Data_Patient, 'None','None')

P_A=data.frame('Patient'=Data_scanner$Patient[which(Data_scanner$Scanner=='A')])


performance_metric=c('AUCpr')


Disccretisation_method=c('FBW','FBN','Med_FBW','Med_FBN','Med_FBW_FBN')

interpolation=c(TRUE)


clinical_F=c('Yes','No')


Type_radiomics=c('OR_Radiomics','TLR_Radiomics')


FS=c('Ranking_FS')


count=0

Test_scanner=c('EF')
Classifier=c('RF','SVM','LR','NB')

for (i_pm in 1:length(performance_metric)) {
  
  for (i_dm in 1:length(Disccretisation_method)) {
    for (i_i in 1:length(interpolation)) {
      
      for (i_tR in 1:length(Type_radiomics)) {
        
        if (Disccretisation_method[i_dm]=='FBW'){
          
          Disccretisation_NR=c(0.05,0.1,0.2,0.5,'ALL')
          
        }
        else{
          Disccretisation_NR=c(32,64,'ALL')
        }
        
        for (i_dmNR in 1:length(Disccretisation_NR)) {
          for (i_cl in 1:length(clinical_F)) {
            for (i_fs in 1:length(FS)) {
              
              if (Disccretisation_method[i_dm]=='FBW' || Disccretisation_method[i_dm]=='FBN'){
              Data_Radiomics=Load_data(General_path,Type_radiomics[i_tR],Data_Patient, Disccretisation_method[i_dm], interpolation[i_i])
              Data_Radiomics=chose_descretisation_nr(Disccretisation_NR[i_dmNR],Data_Radiomics)
              } else{
                if (Disccretisation_method[i_dm]=='Med_FBW'){
                  Data_Radiomics=Load_data(General_path,Type_radiomics[i_tR],Data_Patient, 'FBW', interpolation[i_i])
                  df_r=Data_Radiomics['Patient']
                  Data_Radiomics1=chose_descretisation_nr(0.05,Data_Radiomics)
                  Data_Radiomics1$Patient=NULL
                  Data_Radiomics2=chose_descretisation_nr(0.1,Data_Radiomics)
                  Data_Radiomics2$Patient=NULL
                  Data_Radiomics3=chose_descretisation_nr(0.2,Data_Radiomics)
                  Data_Radiomics3$Patient=NULL
                  Data_Radiomics4=chose_descretisation_nr(0.5,Data_Radiomics)
                  Data_Radiomics4$Patient=NULL
                  col_names=colnames(Data_Radiomics1)
                  for (ind in 1:(dim(Data_Radiomics1)[2])){
                    col_name_f=substring(col_names[ind], 1, nchar(col_names[ind])-5)
                  x_f=apply(data.frame(Data_Radiomics1[,ind],Data_Radiomics2[,ind],Data_Radiomics3[,ind],Data_Radiomics4[,ind]),1,median)
                  df_r[col_name_f]=x_f
                  }
                  Data_Radiomics=df_r
                  
                }
                
                if (Disccretisation_method[i_dm]=='Med_FBN'){
                  Data_Radiomics=Load_data(General_path,Type_radiomics[i_tR],Data_Patient, 'FBN', interpolation[i_i])
                  df_r=Data_Radiomics['Patient']
                  Data_Radiomics1=chose_descretisation_nr(32,Data_Radiomics)
                  Data_Radiomics1$Patient=NULL
                  Data_Radiomics2=chose_descretisation_nr(64,Data_Radiomics)
                  Data_Radiomics2$Patient=NULL
                  
                  col_names=colnames(Data_Radiomics1)
                  for (ind in 1:(dim(Data_Radiomics1)[2])){
                    col_name_f=substring(col_names[ind], 1, nchar(col_names[ind])-5)
                    x_f=apply(data.frame(Data_Radiomics1[,ind],Data_Radiomics2[,ind]),1,median)
                    df_r[col_name_f]=x_f
                  }
                  Data_Radiomics=df_r
                  
                }
                
                if (Disccretisation_method[i_dm]=='Med_FBW_FBN'){
                  Data_Radiomics_FBW=Load_data(General_path,Type_radiomics[i_tR],Data_Patient, 'FBW', interpolation[i_i])
                  Data_Radiomics_FBN=Load_data(General_path,Type_radiomics[i_tR],Data_Patient, 'FBW', interpolation[i_i])
                  
                  df_r=Data_Radiomics_FBW['Patient']
                  Data_Radiomics1=chose_descretisation_nr(0.05,Data_Radiomics_FBW)
                  Data_Radiomics1$Patient=NULL
                  Data_Radiomics2=chose_descretisation_nr(0.1,Data_Radiomics_FBW)
                  Data_Radiomics2$Patient=NULL
                  Data_Radiomics3=chose_descretisation_nr(0.2,Data_Radiomics_FBW)
                  Data_Radiomics3$Patient=NULL
                  Data_Radiomics4=chose_descretisation_nr(0.5,Data_Radiomics_FBW)
                  Data_Radiomics4$Patient=NULL
                  Data_Radiomics5=chose_descretisation_nr(32,Data_Radiomics_FBN)
                  Data_Radiomics5$Patient=NULL
                  Data_Radiomics6=chose_descretisation_nr(64,Data_Radiomics_FBN)
                  Data_Radiomics6$Patient=NULL
                  col_names=colnames(Data_Radiomics1)
                  for (ind in 1:(dim(Data_Radiomics1)[2])){
                    col_name_f=substring(col_names[ind], 1, nchar(col_names[ind])-5)
                    x_f=apply(data.frame(Data_Radiomics1[,ind],Data_Radiomics2[,ind],Data_Radiomics3[,ind],Data_Radiomics4[,ind],Data_Radiomics5[,ind],Data_Radiomics6[,ind]),1,median)
                    df_r[col_name_f]=x_f
                  }
                  Data_Radiomics=df_r
                  
                }
              }
              
              
              
              Data =merge( DATA_C,DATA_output,by="Patient")
              Data =merge( Data,Data_Radiomics,by="Patient")
              
              D_T=Data
              D_T$Time=NULL
              D_T['Outcome']=Data$DFS
              D_T['Patient']=Data$Patient
              
              D_T['Outcome']=as.factor(D_T$Outcome)
              
              
              
              D_T$DFS=NULL
              D_T2=D_T
              D_T2 =merge( D_T2,P_A,by="Patient")
              Out_DT2=D_T2$Outcome
              D_T2$Outcome=NULL
              D_T2$Patient=NULL
              
              D_T2=D_T2[is.finite(colSums(D_T2)) ]
              D_T2 =D_T2[, sapply(data.frame(sapply(D_T2, as.numeric)), var) != 0]
              
              D_T2_IN=D_T2
              
              D_T2_IN['Outcome']=NULL
              
              if(clinical_F[i_cl]=='No'){
                D_T2_IN=D_T2_IN[5:dim(D_T2_IN)[2]]
                
              }
              
        
                TR_Features_MRMR_Forward_pearson=FeatureSelection_MLMR(D_T2_IN,Out_DT2,length(colnames(D_T2_IN)),'up','pearson')
                TR_Features_MRMR_Forward_pearson=colnames(TR_Features_MRMR_Forward_pearson)
                TR_Features_MRMR_Forward_pearson_DF=data.frame('features'=TR_Features_MRMR_Forward_pearson, 'Score1'=order(seq(1, length(colnames(D_T2_IN)), by=1), decreasing = FALSE))
                
                TR_Features_MRMR_Forward_spearman=FeatureSelection_MLMR(D_T2_IN,Out_DT2,length(colnames(D_T2_IN)),'up','spearman')
                TR_Features_MRMR_Forward_spearman=colnames(TR_Features_MRMR_Forward_spearman)
                TR_Features_MRMR_Forward_spearman_DF=data.frame('features'=TR_Features_MRMR_Forward_spearman, 'Score2'=order(seq(1, length(colnames(D_T2_IN)), by=1), decreasing = FALSE))
                print('TR_Features_MRMR_Forward_spearman')
                
                lid <- matrix(0, length(colnames(D_T2_IN)), 1)
                old_f=colnames(D_T2_IN)
                TR_F=MRMR(D_T2_IN, Out_DT2, k = length(colnames(D_T2_IN)))
                for (nn in 1:(length(TR_F$selection))){
                  lid[nn]=TR_F$selection[[nn]]
                }
                TR_Features_MRMR_MI=old_f[lid]
                TR_Features_MRMR_MI_DF=data.frame('features'=TR_Features_MRMR_MI, 'Score5'=order(seq(1, length(colnames(D_T2_IN)), by=1), decreasing = TRUE))
                
                TR_Features_RF_Gini=FeatureSelection_RandomForest(D_T2_IN,1,Out_DT2,'MeanDecreaseGini',length(colnames(D_T2_IN)))
                TR_Features_RF_Gini=colnames(TR_Features_RF_Gini)
                TR_Features_RF_Gini_DF=data.frame('features'=TR_Features_RF_Gini, 'Score6'=order(seq(1, length(colnames(D_T2_IN)), by=1), decreasing = TRUE))
                
                TR_Features_RF_Accuracy=FeatureSelection_RandomForest(D_T2_IN,1,Out_DT2,'MeanDecreaseAccuracy',length(colnames(D_T2_IN)))
                TR_Features_RF_Accuracy=colnames(TR_Features_RF_Accuracy)
                TR_Features_RF_Accuracy_DF=data.frame('features'=TR_Features_RF_Accuracy, 'Score7'=order(seq(1, length(colnames(D_T2_IN)), by=1), decreasing = TRUE))
                print('TR_Features_RF_Accuracy')
                
                
                Merge_Features1=merge(TR_Features_MRMR_Forward_pearson_DF,TR_Features_MRMR_Forward_spearman_DF,by='features')
                Merge_Features4=merge(Merge_Features1,TR_Features_MRMR_MI_DF,by='features')
                Merge_Features5=merge(Merge_Features4,TR_Features_RF_Gini_DF,by='features')
                Merge_Features6=merge(Merge_Features5,TR_Features_RF_Accuracy_DF,by='features')
                Merge_Features_final=Merge_Features6
                Merge_Features_final['Score']=Merge_Features_final['Score1']+Merge_Features_final['Score2']+
                  +Merge_Features_final['Score5']+
                  Merge_Features_final['Score6']+Merge_Features_final['Score7']
                
                Merge_Features_final=Merge_Features_final[order(-Merge_Features_final$Score),]
                
                Features_10=Merge_Features_final['features'][1:10,1]
                Features_5=Merge_Features_final['features'][1:5,1]
              
              name_m_i_pm=paste(performance_metric[i_pm],'_', sep = "_", collapse = NULL)
              name_m_i_fs=paste(FS[i_fs],'_', sep = "_", collapse = NULL)
              name_m_i_dm=paste(Disccretisation_method[i_dm],'_', sep = "_", collapse = NULL)
              name_m_i_i=paste('interpolation_',interpolation[i_i],'_', sep = "_", collapse = NULL)
              name_m_i_cl=paste('Clinical_',clinical_F[i_cl],'_', sep = "_", collapse = NULL)
              name_m_i_tR=paste(Type_radiomics[i_tR],'_', sep = "_", collapse = NULL)
              name_m_idmNR=paste(Disccretisation_NR[i_dmNR],'_',"Combine_FS_","Results_compare_FS.RData", sep = "_", collapse = NULL)
              name_file_Results=paste(name_m_i_pm,name_m_i_fs,name_m_i_dm,name_m_i_i,name_m_i_cl,name_m_i_tR,name_m_idmNR, sep = "_", collapse = NULL)
              print(name_file_Results)
              
              
              Data_DISR=Data
              
              
              Data_DISR['Outcome']=D_T$Outcome
              Data_DISR['Patient']=D_T$Patient
              
              Results_10=sapply(Test_scanner,Performance_Total,Classifier,Data_DISR,Data_scanner,Features_10)
              
              Results_Val_10=Results_10[1,]
              Results_Test_EF_10=Results_10[2,]
              
              RF_Val_10=c(Results_Val_10[[1]][[1]])
              SVM_Val_10=c(Results_Val_10[[1]][[2]])
              LR_Val_10=c(Results_Val_10[[1]][[3]])
              NB_Val_10=c(Results_Val_10[[1]][[4]])
              
              median_Val_10=median(c(Results_Val_10[[1]][[1]],Results_Val_10[[1]][[2]],Results_Val_10[[1]][[3]],Results_Val_10[[1]][[4]]))
              min_Val_10=min(c(Results_Val_10[[1]][[1]],Results_Val_10[[1]][[2]],Results_Val_10[[1]][[3]],Results_Val_10[[1]][[4]]))
              
              RF_test_EF_10=c(Results_Test_EF_10[[1]][[1]])
              SVM_test_EF_10=c(Results_Test_EF_10[[1]][[2]])
              LR_test_EF_10=c(Results_Test_EF_10[[1]][[3]])
              NB_test_EF_10=c(Results_Test_EF_10[[1]][[4]])
              
              median_test_EF_10=median(c(Results_Test_EF_10[[1]][[1]],Results_Test_EF_10[[1]][[2]],Results_Test_EF_10[[1]][[3]],Results_Test_EF_10[[1]][[4]]))
              min_test_EF_10=min(c(Results_Test_EF_10[[1]][[1]],Results_Test_EF_10[[1]][[2]],Results_Test_EF_10[[1]][[3]],Results_Test_EF_10[[1]][[4]]))
              
              Results_10=sapply(c('H'),Performance_Total,Classifier,Data_DISR,Data_scanner,Features_10)
              Results_Test_H_10=Results_10[2,]
              
              RF_test_H_10=c(Results_Test_H_10[[1]][[1]])
              SVM_test_H_10=c(Results_Test_H_10[[1]][[2]])
              LR_test_H_10=c(Results_Test_H_10[[1]][[3]])
              NB_test_H_10=c(Results_Test_H_10[[1]][[4]])
              
              
              median_test_H_10=median(c(Results_Test_H_10[[1]][[1]],Results_Test_H_10[[1]][[2]],Results_Test_H_10[[1]][[3]],Results_Test_H_10[[1]][[4]]))
              min_test_H_10=min(c(Results_Test_H_10[[1]][[1]],Results_Test_H_10[[1]][[2]],Results_Test_H_10[[1]][[3]],Results_Test_H_10[[1]][[4]]))
              
              
              Results_5=sapply(Test_scanner,Performance_Total,Classifier,Data_DISR,Data_scanner,Features_5)
              
              Results_Val_5=Results_5[1,]
              Results_Test_EF_5=Results_5[2,]
              
              RF_Val_5=c(Results_Val_5[[1]][[1]])
              SVM_Val_5=c(Results_Val_5[[1]][[2]])
              LR_Val_5=c(Results_Val_5[[1]][[3]])
              NB_Val_5=c(Results_Val_5[[1]][[4]])
              
              median_Val_5=median(c(Results_Val_5[[1]][[1]],Results_Val_5[[1]][[2]],Results_Val_5[[1]][[3]],Results_Val_5[[1]][[4]]))
              min_Val_5=min(c(Results_Val_5[[1]][[1]],Results_Val_5[[1]][[2]],Results_Val_5[[1]][[3]],Results_Val_5[[1]][[4]]))
              
              RF_test_EF_5=c(Results_Test_EF_5[[1]][[1]])
              SVM_test_EF_5=c(Results_Test_EF_5[[1]][[2]])
              LR_test_EF_5=c(Results_Test_EF_5[[1]][[3]])
              NB_test_EF_5=c(Results_Test_EF_5[[1]][[4]])
              
              median_test_EF_5=median(c(Results_Test_EF_5[[1]][[1]],Results_Test_EF_5[[1]][[2]],Results_Test_EF_5[[1]][[3]],Results_Test_EF_5[[1]][[4]]))
              min_test_EF_5=min(c(Results_Test_EF_5[[1]][[1]],Results_Test_EF_5[[1]][[2]],Results_Test_EF_5[[1]][[3]],Results_Test_EF_5[[1]][[4]]))
              
              Results_5=sapply(c('H'),Performance_Total,Classifier,Data_DISR,Data_scanner,Features_5)
              Results_Test_H_5=Results_5[2,]
              
              RF_test_H_5=c(Results_Test_H_5[[1]][[1]])
              SVM_test_H_5=c(Results_Test_H_5[[1]][[2]])
              LR_test_H_5=c(Results_Test_H_5[[1]][[3]])
              NB_test_H_5=c(Results_Test_H_5[[1]][[4]])
              
              median_test_H_5=median(c(Results_Test_H_5[[1]][[1]],Results_Test_H_5[[1]][[2]],Results_Test_H_5[[1]][[3]],Results_Test_H_5[[1]][[4]]))
              min_test_H_5=min(c(Results_Test_H_5[[1]][[1]],Results_Test_H_5[[1]][[2]],Results_Test_H_5[[1]][[3]],Results_Test_H_5[[1]][[4]]))
              
              if (count==0){
                Results_df_10=data.frame('Name'=name_file_Results, 'Features'=paste(Features_10, collapse = ';'),'median_val'=median_Val_10,'min_val'=min_Val_10, 'RF_Val'=RF_Val_10,'SVM_Val'=SVM_Val_10, 'LR_Val'=LR_Val_10,'NB_Val'=NB_Val_10, 'median_test_EF'=median_test_EF_10, 'min_test_EF'=min_test_EF_10,'RF_test_EF'=RF_test_EF_10,'SVM_test_EF'=SVM_test_EF_10, 'LR_test_EF'=LR_test_EF_10,'NB_test_EF'=NB_test_EF_10, 'median_test_H'=median_test_H_10, 'min_test_H'=min_test_H_10, 'RF_test_H'=RF_test_H_10,'SVM_test_H'=SVM_test_H_10, 'LR_test_H'=LR_test_H_10,'NB_test_H'=NB_test_H_10)
                Results_df_5=data.frame('Name'=name_file_Results, 'Features'=paste(Features_5, collapse = ';'),'median_val'=median_Val_5,'min_val'=min_Val_5, 'RF_Val'=RF_Val_5,'SVM_Val'=SVM_Val_5, 'LR_Val'=LR_Val_5,'NB_Val'=NB_Val_5, 'median_test_EF'=median_test_EF_5, 'min_test_EF'=min_test_EF_5,'RF_test_EF'=RF_test_EF_5,'SVM_test_EF'=SVM_test_EF_5, 'LR_test_EF'=LR_test_EF_5,'NB_test_EF'=NB_test_EF_5, 'median_test_H'=median_test_H_5, 'min_test_H'=min_test_H_5, 'RF_test_H'=RF_test_H_5,'SVM_test_H'=SVM_test_H_5, 'LR_test_H'=LR_test_H_5,'NB_test_H'=NB_test_H_5)
                
              }else{
                Results_df_10_new=data.frame('Name'=name_file_Results, 'Features'=paste(Features_10, collapse = ';'),'median_val'=median_Val_10,'min_val'=min_Val_10, 'RF_Val'=RF_Val_10,'SVM_Val'=SVM_Val_10, 'LR_Val'=LR_Val_10,'NB_Val'=NB_Val_10, 'median_test_EF'=median_test_EF_10, 'min_test_EF'=min_test_EF_10,'RF_test_EF'=RF_test_EF_10,'SVM_test_EF'=SVM_test_EF_10, 'LR_test_EF'=LR_test_EF_10,'NB_test_EF'=NB_test_EF_10, 'median_test_H'=median_test_H_10, 'min_test_H'=min_test_H_10, 'RF_test_H'=RF_test_H_10,'SVM_test_H'=SVM_test_H_10, 'LR_test_H'=LR_test_H_10,'NB_test_H'=NB_test_H_10)
                Results_df_10=rbind(Results_df_10,Results_df_10_new)
                
                Results_df_5_new=data.frame('Name'=name_file_Results, 'Features'=paste(Features_5, collapse = ';'),'median_val'=median_Val_5,'min_val'=min_Val_5, 'RF_Val'=RF_Val_5,'SVM_Val'=SVM_Val_5, 'LR_Val'=LR_Val_5,'NB_Val'=NB_Val_5, 'median_test_EF'=median_test_EF_5, 'min_test_EF'=min_test_EF_5,'RF_test_EF'=RF_test_EF_5,'SVM_test_EF'=SVM_test_EF_5, 'LR_test_EF'=LR_test_EF_5,'NB_test_EF'=NB_test_EF_5, 'median_test_H'=median_test_H_5, 'min_test_H'=min_test_H_5, 'RF_test_H'=RF_test_H_5,'SVM_test_H'=SVM_test_H_5, 'LR_test_H'=LR_test_H_5,'NB_test_H'=NB_test_H_5)
                
                Results_df_5=rbind(Results_df_5,Results_df_5_new)
              }
              
              count=count+1
              
            }
            
            
            
          }
        }
      }
    }
  }
}

