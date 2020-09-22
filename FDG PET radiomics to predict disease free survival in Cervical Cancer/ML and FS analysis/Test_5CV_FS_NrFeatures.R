Test_5CV_FS_NrFeatures= function (TestSet, features,Set, folding, folds, Var3, Z){

#####################################
# Test with k fold cross validation

TestSet_out= TestSet$Outcome
TestSet2=TestSet

TestSet2= TestSet2[c(features)]

TestSet2['Outcome']=TestSet_out

AUC_TESTE <- matrix(0, folds, 1)
fscore_TESTE <- matrix(0, folds, 1)
prauc_TESTE <- matrix(0, folds, 1)
precision_TESTE <- matrix(0, folds, 1)
recall_TESTE <- matrix(0, folds, 1)
Acc <- matrix(0, folds, 1)
t_cross <- matrix(0, folds, 1)
sensitivity_TE<- matrix(0, folds, 1)
sensitivity_TE2<- matrix(0, folds, 1)
specificity_TE2<- matrix(0, folds, 1)

# For each k fold
for (i in 1:folds){
  
  
  # Define train and validation fold
  vTR=list(folding$train$'1'$idx, folding$train$'2'$idx, folding$train$'3'$idx, folding$train$'4'$idx,
           folding$train$'5'$idx)
  vTE=list(folding$test$'1'$idx, folding$test$'2'$idx, folding$test$'3'$idx, folding$test$'4'$idx,
           folding$test$'5'$idx)
  
  TR=Set[vTR[i][[1]],]
  TE=Set[vTE[i][[1]],]
  
  
  TR_Out=TR$Outcome
 
  TR=TR[c(features)]
 
  TR['Outcome']=TR_Out
  
  TE_out= TE$Outcome
  
  TE= TE[c(features)]
  
  TE['Outcome']=TE_out
  
  
  # Train and test with RF
  set.seed(Z)
  if (Var3==1){
    model<- randomForest(Outcome ~ ., data = TR, importance = TRUE)
   
    predTest <- predict(model, TestSet2, type = "prob")
    predTest2=prediction(predTest[,2], TestSet2$Outcome)
  }
  if (Var3==2){
    model<- svm(Outcome ~ ., data = TR, Kernel=kernel_svm,probability = TRUE)
   
    predTest <- predict(model, TestSet2, type = "prob",probability = TRUE)
    predTest= attr(predTest, "probabilities")
    predTest2=prediction(predTest[,2], TestSet2$Outcome)
  }
  if (Var3==3){
    model <- glm(as.formula(paste('Outcome', paste(features, collapse=" + "), sep=" ~ ")), 
                 data = TR, 
                 family = binomial)
   
    predTest <- predict(model, TestSet2, type = "response")
    predTest2=prediction(predTest, TestSet2$Outcome)
  }
  if (Var3==4){
    model<- naiveBayes(Outcome ~ ., data = TR)
   
    predTest <- predict(model, TestSet2, type = "raw")
    predTest2=prediction(predTest[,2], TestSet2$Outcome)
  }
  
 
  
  ROC_TEST=performance(predTest2,"sens", "spec")
  
  AUC_TESTE[i]=performance(predTest2,measure = "auc")@y.values[[1]]
  
  YI=ROC_TEST@y.values[[1]]+ROC_TEST@x.values[[1]]-1
  ind_YI= match(max(YI),YI)
  THRESH= ROC_TEST@alpha.values[[1]][ind_YI]
 
  

  
  sensitivity_TE[i]=ROC_TEST@y.values[[1]][ind_YI]
  
  
  
  ROC_TEST_pr=performance(predTest2,"prec", "rec")
  
  
  
  YI_pr=ROC_TEST_pr@y.values[[1]]+ROC_TEST_pr@x.values[[1]]-1
  ind_YI_pr= match(max(YI_pr[2:length(YI_pr)]),YI_pr)
  THRESH_pr= ROC_TEST_pr@alpha.values[[1]][ind_YI_pr]
  
  a=(performance(predTest2,"prec", "rec")@y.values[[1]])
  a=a[2:length(a)]
  b=diff(performance(predTest2,"prec", "rec")@x.values[[1]])
  prauc=sum(a*b)
  
  ind_YI_pr=length(which(as.numeric(predTest2@cutoffs[[1]])>0.5))
  
  
  tp=predTest2@tp[[1]][ind_YI_pr]
  tn=predTest2@tn[[1]][ind_YI_pr]
  fn=predTest2@fn[[1]][ind_YI_pr]
  fp=predTest2@fp[[1]][ind_YI_pr]
  
  
  beta=1
  
  if (tp==0){
    recall=0
    precision=0
    fscore=0
  }else{
    recall=tp/(tp+fn)
    precision=tp/(tp+fp)
  fscore=((1+(beta^2))*(precision*recall))/(((beta^2)*precision)+recall)
  }
 
  
  fscore_TESTE[i]=fscore
  prauc_TESTE[i]=prauc
  precision_TESTE[i]=precision
  recall_TESTE[i]=recall
 
  
  
  if (Var3==3)
  {
    predTest_THRESH=predTest
  }else{
    predTest_THRESH=predTest[,2]
    
  }
  predTest_THRESH[which(predTest<THRESH_pr)]=1#THRESH
  predTest_THRESH[which(predTest>=THRESH_pr)]=0
  Acc[i]=Accuracy(predTest_THRESH, TestSet2$Outcome)
  
  
  TP=length(which(((predTest_THRESH==1)&(TestSet2$Outcome==1))==TRUE))
  TN=length(which(((predTest_THRESH==0)&(TestSet2$Outcome==0))==TRUE))
  FP=length(which(((predTest_THRESH==1)&(TestSet2$Outcome==0))==TRUE))
  FN=length(which(((predTest_THRESH==0)&(TestSet2$Outcome==1))==TRUE))
  
  
  sensitivity_TE2[i]=TP/(TP+FN)
  specificity_TE2[i]=TN/(TN+FP)
  
  
}

CI_down_TEST=min(AUC_TESTE)
CI_up_TEST=max(AUC_TESTE)
mean_TEST=mean(AUC_TESTE)

CI_down_TEST_fscore=min(fscore_TESTE)
CI_up_TEST_fscore=max(fscore_TESTE)
mean_TEST_fscore=mean(fscore_TESTE)
mean_TEST_prauc=mean(prauc_TESTE)
mean_TEST_precision=mean(precision_TESTE)
mean_TEST_recall=mean(recall_TESTE)

CI_down_sensitivity_TEST=min(sensitivity_TE)
CI_up_sensitivity_TEST=max(sensitivity_TE)
mean_sensitivity_TEST=mean(sensitivity_TE)

return(c(CI_down_TEST,CI_up_TEST,mean_TEST, CI_down_TEST_fscore,CI_up_TEST_fscore,mean_TEST_fscore, mean_TEST_precision,mean_TEST_recall,mean_TEST_prauc))

}