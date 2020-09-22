CV_Select_FS_NrFeatures= function (Z,Set,N, folds, FS, Var3, folding, criteria,cri,v_f,beta){
 
  
  sensitivity_TRAIN <- matrix(0, dim(Data)[2], folds)
  sensitivity_CROSS <- matrix(0, dim(Data)[2], folds)
  
  specificity_TRAIN <- matrix(0, dim(Data)[2], folds)
  specificity_CROSS <- matrix(0, dim(Data)[2], folds)
  
  AUC_Cross_R_F_cross <- matrix(0, dim(Data)[2], folds)
  AUC_TRAIN_R_F_cross <- matrix(0, dim(Data)[2], folds)
  
  fscore_Cross_R_F_cross <- matrix(0, dim(Data)[2], folds)
  fscore_TRAIN_R_F_cross <- matrix(0, dim(Data)[2], folds)
  
  prauc_Cross_R_F_cross <- matrix(0, dim(Data)[2], folds)
  
  Balance_TE<- matrix(0, dim(Data)[2], folds)
  Balance_TR<- matrix(0, dim(Data)[2], folds)
  
  # Test para numeros diferentes de features e usar 5 fold cross validation to test these features
  sensitivity_TRAIN_R_F <- matrix(0, N, 1)
  sensitivity_CROSS_R_F <- matrix(0, N, 1)
  
  specificity_TRAIN_R_F <- matrix(0, N, 1)
  specificity_CROSS_R_F <- matrix(0, N, 1)
  
  CI_down_sensitivity <- matrix(0, N, 1)
  CI_up_sensitivity <- matrix(0, N, 1)
  CI_down_specificity <- matrix(0, N, 1)
  CI_up_specificity <- matrix(0, N, 1)
  
  AUC_TRAIN_R_F <- matrix(0, N, 1)
  AUC_Cross_R_F <- matrix(0, N, 1)
  
  fscore_TRAIN_R_F <- matrix(0, N, 1)
  fscore_Cross_R_F <- matrix(0, N, 1)
  
  prauc_Cross_R_F <- matrix(0, N, 1)
  
  CI_down <- matrix(0, N, 1)
  CI_up <- matrix(0, N, 1)
  
  CI_down_fscore <- matrix(0, N, 1)
  CI_up_fscore <- matrix(0, N, 1)
  
  # for each number of features
  for (j in 1:N)
  {
    
    NR_Features=v_f[j]
    
    TR_Features=Set
    
    TR_Features['Outcome']=NULL
    
    
    
    #Select features according to method
    if ((FS=='MRMR_pearson')|(FS=='MRMR_spearman')){
      
      TR_Features=FeatureSelection_MLMR(TR_Features,Set$Outcome,NR_Features,direction,FS)
      features=colnames(TR_Features)
      
    }
    if (FS=='RF'){
      TR_Features=FeatureSelection_RandomForest(TR_Features,Z,Set$Outcome,criteria,NR_Features)
      features=colnames(TR_Features)
    }
    
    lid <- matrix(0, NR_Features, 1) 
    
    if (FS=='MRMR_MI'){
      old_f=colnames(TR_Features)
      TR_F=MRMR(TR_Features, Set$Outcome, k = NR_Features)
      for (nn in 1:(length(TR_F$selection))){
        lid[nn]=TR_F$selection[[nn]]
      }
      features=old_f[lid]
    }
    
    
    # For each k fold
    for (i in 1:folds)
    {
      
   
      # Define train and validation fold
      vTR=list(folding$train$'1'$idx, folding$train$'2'$idx, folding$train$'3'$idx, folding$train$'4'$idx,
               folding$train$'5'$idx)
      vTE=list(folding$test$'1'$idx, folding$test$'2'$idx, folding$test$'3'$idx, folding$test$'4'$idx,
               folding$test$'5'$idx)
      
      
      TR=Set[vTR[i][[1]],]
      TE=Set[vTE[i][[1]],]
      
      Outcome_TR=TR['Outcome']
      Outcome_TE=TE['Outcome']
      
     
      TR=TR[c(features)]
      TE= TE[c(features)]
    
      TR['Outcome']=Outcome_TR
      TE['Outcome']=Outcome_TE
     
      
      TR =TR[, sapply(data.frame(sapply(TR, as.numeric)), var) != 0]
     
      
      if (typeof(TR)!="integer"){
        features=colnames(TR)[(which(colnames(TR)!="Outcome"))]
      # Train and test 
      set.seed(Z)
      if (Var3==1){
        model<- randomForest(Outcome ~ ., data = TR, importance = TRUE)
        predTrain <- predict(model, TR, type = "prob")
        predTrain2 <- prediction(predTrain[,2], TR$Outcome)
        predCross <- predict(model, TE, type = "prob")
        predCross2 <- prediction(predCross[,2], TE$Outcome)
      }
      if (Var3==2){
        model<- svm(Outcome ~ ., data = TR, Kernel=kernel_svm,probability = TRUE)
        predTrain <- predict(model, TR, type = "prob",probability = TRUE)
        predTrain= attr(predTrain, "probabilities")
        predTrain2 <- prediction(predTrain[,2], TR$Outcome)
        predCross <- predict(model, TE, type = "prob",probability = TRUE)
        predCross= attr(predCross, "probabilities")
        predCross2 <- prediction(predCross[,2], TE$Outcome)
      }
      if (Var3==3){
        model <- glm(as.formula(paste('Outcome', paste(features, collapse=" + "), sep=" ~ ")), 
                     data = TR, 
                     family = binomial)
        predTrain <- predict(model, TR, type = "response")
        predTrain2=prediction(predTrain, TR$Outcome)
        predCross <- predict(model, TE, type = "response")
        predCross2=prediction(predCross, TE$Outcome)
      }
      if (Var3==4){
        model<- naiveBayes(Outcome ~ ., data = TR)
        predTrain <- predict(model, TR, type = "raw")
        predTrain2 <- prediction(predTrain[,2], TR$Outcome)
        predCross <- predict(model, TE, type = "raw")
        predCross2 <- prediction(predCross[,2], TE$Outcome)
        
      }
      }else{
       
        next
      }
      
      
      AUC_TRAIN=performance(predTrain2,measure = "auc")@y.values[[1]]
      ROC_TRAIN=performance(predTrain2,"sens", "spec")
      
      
      YI=ROC_TRAIN@y.values[[1]]+ROC_TRAIN@x.values[[1]]-1
      ind_YI= match(max(YI),YI)
      THRESH= ROC_TRAIN@alpha.values[[1]][ind_YI]
      
      
      sensitivity_TRAIN[j,i]=ROC_TRAIN@y.values[[1]][ind_YI]
      specificity_TRAIN[j,i]=ROC_TRAIN@x.values[[1]][ind_YI]
      
      ROC_TRAIN_pr=performance(predTrain2,"prec", "rec")
      
      YI_pr=ROC_TRAIN_pr@y.values[[1]]+ROC_TRAIN_pr@x.values[[1]]-1
      ind_YI_pr= match(max(YI_pr[2:length(YI_pr)]),YI_pr)
      THRESH_pr= ROC_TRAIN_pr@alpha.values[[1]][ind_YI_pr]
      
      ind_YI_pr=length(which(as.numeric(predTrain2@cutoffs[[1]])>0.5))
      
      tp=predTrain2@tp[[1]][ind_YI_pr]
      tn=predTrain2@tn[[1]][ind_YI_pr]
      fn=predTrain2@fn[[1]][ind_YI_pr]
      fp=predTrain2@fp[[1]][ind_YI_pr]
      
      
      
      
      if (tp==0){
        fscore_TRAIN=0
        recall=0
        precision=0
      }else{
        recall=tp/(tp+fn)
        precision=tp/(tp+fp)
      fscore_TRAIN=((1+(beta^2))*(precision*recall))/(((beta^2)*precision)+recall)
      }
      
      if (class(try(performance( predCross2,measure = "auc")@y.values[[1]], silent=T))=="try-error"){
        AUC_Cross=0
        sensitivity_CROSS[j,i]=0
        specificity_CROSS[j,i]=0
        fscore_Cross=0
        prauc_Cross=0
        
        }
      else{
      AUC_Cross =performance( predCross2,measure = "auc")@y.values[[1]]
      ROC_Cross=performance( predCross2,"sens", "spec")
      
      YI=ROC_Cross@y.values[[1]]+ROC_Cross@x.values[[1]]-1
      ind_YI= match(max(YI),YI)
      THRESH= ROC_Cross@alpha.values[[1]][ind_YI]
      
      sensitivity_CROSS[j,i]=ROC_Cross@y.values[[1]][ind_YI]
      specificity_CROSS[j,i]=ROC_Cross@x.values[[1]][ind_YI]
      
      ROC_Cross_pr=performance(predCross2,"prec", "rec")
      
      YI_pr=ROC_Cross_pr@y.values[[1]]+ROC_Cross_pr@x.values[[1]]-1
      ind_YI_pr= match(max(YI_pr[2:length(YI_pr)]),YI_pr)
      THRESH_pr= ROC_Cross_pr@alpha.values[[1]][ind_YI_pr]
      
      ind_YI_pr=length(which(as.numeric(predCross2@cutoffs[[1]])>0.5))
      
      tp=predCross2@tp[[1]][ind_YI_pr]
      tn=predCross2@tn[[1]][ind_YI_pr]
      fn=predCross2@fn[[1]][ind_YI_pr]
      fp=predCross2@fp[[1]][ind_YI_pr]
      
     
      
      
      if (tp==0){
        fscore_Cross=0
        recall=0
        precision=0
      }else{
        recall=tp/(tp+fn)
        precision=tp/(tp+fp)
      fscore_Cross=((1+(beta^2))*(precision*recall))/(((beta^2)*precision)+recall)
      }
      
      a=(performance(predCross2,"prec", "rec")@y.values[[1]])
      a=a[2:length(a)]
      b=diff(performance(predCross2,"prec", "rec")@x.values[[1]])
      prauc_Cross=sum(a*b)
      
      
      AUC_Cross_R_F_cross[j,i]=AUC_Cross
      AUC_TRAIN_R_F_cross[j,i]=AUC_TRAIN
      
      fscore_Cross_R_F_cross[j,i]=fscore_Cross
      fscore_TRAIN_R_F_cross[j,i]=fscore_TRAIN
      
      prauc_Cross_R_F_cross[j,i]=prauc_Cross
      
      }
      
      
    }
    
    #AUC for each number of features used 
    sensitivity_TRAIN_R_F[j]= mean(sensitivity_TRAIN[j,1:folds])
    sensitivity_CROSS_R_F[j]= mean(sensitivity_CROSS[j,1:folds])
    specificity_TRAIN_R_F[j]= mean(specificity_TRAIN[j,1:folds])
    specificity_CROSS_R_F[j]= mean(specificity_CROSS[j,1:folds])
    
    # Validation AUC confidence intervals for each number of features used 
    CI_down_sensitivity[j]=min(sensitivity_CROSS[j,1:folds])
    CI_up_sensitivity[j]=max(sensitivity_CROSS[j,1:folds])
    CI_down_specificity[j]=min(specificity_CROSS[j,1:folds])
    CI_up_specificity[j]=max(specificity_CROSS[j,1:folds])
    
    #AUC for each number of features used 
    AUC_TRAIN_R_F[j]= mean(AUC_TRAIN_R_F_cross[j,1:folds])
    AUC_Cross_R_F[j]= mean(AUC_Cross_R_F_cross[j,1:folds])
    
    #fscore for each number of features used 
    fscore_TRAIN_R_F[j]= mean(fscore_TRAIN_R_F_cross[j,1:folds])
    fscore_Cross_R_F[j]= mean(fscore_Cross_R_F_cross[j,1:folds])
    
    prauc_Cross_R_F[j]= mean(prauc_Cross_R_F_cross[j,1:folds])
    
    # Validation AUC confidence intervals for each number of features used 
    CI_down[j]=min(AUC_Cross_R_F_cross[j,1:folds])
    CI_up[j]=max(AUC_Cross_R_F_cross[j,1:folds])
    
    # Validation AUC confidence intervals for each number of features used 
    CI_down_fscore[j]=min(fscore_Cross_R_F_cross[j,1:folds])
    CI_up_fscore[j]=max(fscore_Cross_R_F_cross[j,1:folds])
    
  }
  
  


CI_DF=data.frame('Features Number'= v_f, 'Mean AUC'= AUC_Cross_R_F, 'CI_down'=CI_down, 'CI_up'=CI_up,'Mean sensitivity'= sensitivity_CROSS_R_F, 'CI_down sen'=CI_down_sensitivity, 'CI_up sen'=CI_up_sensitivity, 'Mean specificity'= specificity_CROSS_R_F, 'CI_down spe'=CI_down_specificity, 'CI_up spe'=CI_up_specificity, 'Mean fscore'= fscore_Cross_R_F, 'CI_down_fscore'=CI_down_fscore, 'CI_up_fscore'=CI_up_fscore, 'Mean prauc'= prauc_Cross_R_F )

C1_DF_10=CI_DF[1:10,]
#Order the desired variable
C1_DF_10_auc_order=C1_DF_10[order(C1_DF_10[cri]),]


NR_Features=C1_DF_10_auc_order$Features.Number[dim(C1_DF_10_auc_order)[1]]


return(c(NR_Features,C1_DF_10_auc_order))

}
