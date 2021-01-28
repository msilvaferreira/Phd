Performance_TestData <- function(Test_scanner,Data,Data_scanner){
  
  Z=1
  ratio=0.8
  
  set.seed(Z)
  
  D_A=Data_scanner[which(Data_scanner$Scanner=='A'),]
  
  D_EF=Data_scanner[which(Data_scanner$Scanner=='EF'),]
  
  D_H=Data_scanner[which(Data_scanner$Scanner=='H'),]
  
  set.seed(Z)
  train_A <- sample(nrow(D_A), ratio*nrow(D_A), replace = FALSE)#
  Set_A <- D_A[train_A,]
  TestSet_A <- D_A[-train_A,]
  
  set.seed(Z)
  train_EF <- sample(nrow(D_EF), ratio*nrow(D_EF), replace = FALSE)#
  Set_EF <- D_EF[train_EF,]
  TestSet_EF <- D_EF[-train_EF,]
  
  TestSet_H <- D_H
  
  
  #NEW
  if (Test_scanner=='A'){
    Set=rbind(D_EF)
    Set$Institute=NULL
    Set$Scanner=NULL
  }
  
  if (Test_scanner=='EF'){
    Set=rbind(D_A)
    Set$Institute=NULL
    Set$Scanner=NULL
  }
  if (Test_scanner=='H'){
    Set=rbind(D_A)
    Set$Institute=NULL
    Set$Scanner=NULL
  }
  
  
  Set=merge(Set,Data,by='Patient')
  
  
  Outcome_Set=Set$Outcome
  
  Z=1
  folds=5
  set.seed(Z)
  
  
  
  folding_A <- crossv_kfold(D_A,k=folds)
  vTR_A=list(D_A[folding_A$train$'1'$idx,], D_A[folding_A$train$'2'$idx,], D_A[folding_A$train$'3'$idx,], D_A[folding_A$train$'4'$idx,],
                                                    D_A[folding_A$train$'5'$idx,])
  vTE_A=list(D_A[folding_A$test$'1'$idx,], D_A[folding_A$test$'2'$idx,], D_A[folding_A$test$'3'$idx,], D_A[folding_A$test$'4'$idx,],
             D_A[folding_A$test$'5'$idx,])

  folding_EF <- crossv_kfold(D_EF,k=folds)
  vTR_EF=list(D_EF[folding_EF$train$'1'$idx,], D_EF[folding_EF$train$'2'$idx,], D_EF[folding_EF$train$'3'$idx,], D_EF[folding_EF$train$'4'$idx,],
              D_EF[folding_EF$train$'5'$idx,])
  vTE_EF=list(D_EF[folding_EF$test$'1'$idx,], D_EF[folding_EF$test$'2'$idx,], D_EF[folding_EF$test$'3'$idx,], D_EF[folding_EF$test$'4'$idx,],
              D_EF[folding_EF$test$'5'$idx,])

  
  if (Test_scanner=='EF'){
    TestSet=rbind(Set_EF,TestSet_EF) 
    vTR=list(rbind(vTR_A[[1]]),rbind(vTR_A[[2]]),rbind(vTR_A[[3]]),rbind(vTR_A[[4]]),rbind(vTR_A[[5]]))
    vTE=list((vTE_A[[1]]),(vTE_A[[2]]),(vTE_A[[3]]),(vTE_A[[4]]),(vTE_A[[5]]))
    
  }
  if (Test_scanner=='H'){
    TestSet=TestSet_H
    vTR=list(rbind(vTR_A[[1]]),rbind(vTR_A[[2]]),rbind(vTR_A[[3]]),rbind(vTR_A[[4]]),rbind(vTR_A[[5]]))
    vTE=list((vTE_A[[1]]),(vTE_A[[2]]),(vTE_A[[3]]),(vTE_A[[4]]),(vTE_A[[5]]))
  }
  
  TestSet$Institute=NULL
  TestSet$Scanner=NULL
  
  TestSet=merge(TestSet,Data,by='Patient')
  
  Outcome_TestSet=TestSet$Outcome
  
 
  
  return(list(Set,Outcome_Set,TestSet,Outcome_TestSet,vTR,vTE))
}



