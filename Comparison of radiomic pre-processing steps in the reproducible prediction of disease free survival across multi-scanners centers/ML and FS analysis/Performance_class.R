Performance_class <- function(Classifier,Set, Output_Set, TestSet, Output_TestSet,Features, VTR, VTE,Data ) {
  
  folds=length(VTE)
  
  beta=2
  
  AUC_TEST_val=matrix(0,folds,1)
  prauc_TEST_val=matrix(0,folds,1)
  fscore_TEST_val=matrix(0,folds,1)
  recall_TEST_val=matrix(0,folds,1)
  THRESH_val=matrix(0,folds,1)
  BA_val=matrix(0,folds,1)
  
  #fold CV on train set    
  for (ival in 1:folds){
    
    TR_val=merge(VTR[ival],Data,by='Patient')
    TR_val_Outcome=TR_val$Outcome
    TR_val=TR_val[Features]
    TR_val =TR_val[ sapply(data.frame(sapply(TR_val, as.numeric)), var) != 0]
    l_TRVAL=length(TR_val)
    TR_val['Outcome']=TR_val_Outcome
    
    
    TE_val=merge(VTE[ival],Data,by='Patient')
    TE_val_Outcome=TE_val$Outcome
    TE_val=TE_val[Features]
    TE_val['Outcome']=TE_val_Outcome
    
   
    if (l_TRVAL==length(Features)){
    
    if (Classifier=='RF'){
      model<- randomForest(Outcome ~ ., data = TR_val, importance = TRUE)
      predTest <- predict(model, TE_val, type = "prob")
      predTest2 <- prediction(predTest[,2], TE_val_Outcome)
    }
    if (Classifier=='SVM'){
      model<- svm(Outcome ~ ., data = TR_val, Kernel='radial',probability = TRUE)
      predTest <- predict(model, TE_val, type = "prob",probability = TRUE)
      predTest= attr(predTest, "probabilities")
      predTest2 <- prediction(predTest[,2], TE_val_Outcome)
    }
    if (Classifier=='LR'){
      model <- glm(as.formula(paste('Outcome', paste(Features, collapse=" + "), sep=" ~ ")), 
                   data = TR_val, 
                   family = binomial)
      predTest <- predict(model, TE_val, type = "response")
      predTest2=prediction(predTest, TE_val_Outcome)
    }
    if (Classifier=='NB'){
      model<- naiveBayes(Outcome ~ ., data = TR_val)
      predTest <- predict(model, TE_val, type = "raw")
      predTest2 <- prediction(predTest[,2], TE_val_Outcome)
      
    }
    
    
    
    a=(performance(predTest2,"prec", "rec")@y.values[[1]])
    l_a=(length(a))
    a=a[2:l_a]
    b=diff(performance(predTest2,"prec", "rec")@x.values[[1]])
    prauc_TEST_val[ival]=sum(a*b)
    
    
  }
 
  
 }
  
  Set =Set[ sapply(data.frame(sapply(Set, as.numeric)), var) != 0]
  l_TRVAL=length(Set)
  
  #Train on all train set and test on test Set
  Set['Outcome']=Output_Set
  
  if (l_TRVAL==length(Features)){
  if (Classifier=='RF'){
    model<- randomForest(Outcome ~ ., data = Set, importance = TRUE)
    predTest <- predict(model, TestSet, type = "prob")
    predTest2 <- prediction(predTest[,2], Output_TestSet)
  }
  if (Classifier=='SVM'){
    model<- svm(Outcome ~ ., data = Set, Kernel='radial',probability = TRUE)
    predTest <- predict(model, TestSet, type = "prob",probability = TRUE)
    predTest= attr(predTest, "probabilities")
    predTest2 <- prediction(predTest[,2], Output_TestSet)
  }
  if (Classifier=='LR'){
    model <- glm(as.formula(paste('Outcome', paste(Features, collapse=" + "), sep=" ~ ")), 
                 data = Set, 
                 family = binomial)
    predTest <- predict(model, TestSet, type = "response")
    predTest2=prediction(predTest, Output_TestSet)
  }
  if (Classifier=='NB'){
    model<- naiveBayes(Outcome ~ ., data = Set)
    predTest <- predict(model, TestSet, type = "raw")
    predTest2 <- prediction(predTest[,2], Output_TestSet)
    
  }
  
  
  
  a=(performance(predTest2,"prec", "rec")@y.values[[1]])
  
  l_a=(length(a))
  a=a[2:l_a]
  b=diff(performance(predTest2,"prec", "rec")@x.values[[1]])
  prauc_TEST=sum(a*b)
  
 
   
  
  }
  
  

  return(c(prauc_TEST_val_mean,prauc_TEST))
  
}   