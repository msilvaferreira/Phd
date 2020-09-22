#######################################################
## Train with bootstrap and test with Test set ########
#######################################################

#Function to be used by boot. For each random sample of boot training data is used as input to a classifier to predict DFS
#Input:
#data data frame with training data
#indices=Number of bootstrap repetitions
#TEST= data used to test
#var3 =string containing the classifier
#features= vector with features to be used in classification
#Output:
#AUC_TRAIN
#AUC_TEST
#sensitivity_TRAIN
#sensitivity_TEST
#specificity_TRAIN
#specificity_TEST
#THRESH = Youden Index threshold of the ROC of the training data
#Acc= accuracy
#THRESH_TEST = Youden Index threshold of the ROC of the test data


my_RF_to_BOOTSTRAP_Sensitivity_AUC <- function(data, indices, TEST,Var3,features,beta){
  dt<-data[indices,]
  groupvars=features
  
  set.seed(1)
  if (Var3==1){
    model<- randomForest(Outcome ~ ., data = dt, importance = TRUE)
    predTrain <- predict(model, dt, type = "prob")
    predTrain2=prediction(predTrain[,2], dt$Outcome)
    predTest <- predict(model, TEST, type = "prob")
    predTest2=prediction(predTest[,2], TEST$Outcome)
  }
  if (Var3==2){
    model<- svm(Outcome ~ ., data = dt,Kernel=kernel_svm,probability = TRUE)
    predTrain <- predict(model, dt, type = "prob", probability = TRUE)
    predTrain= attr(predTrain, "probabilities")
    predTrain2=prediction(predTrain[,2], dt$Outcome)
    predTest <- predict(model, TEST, type = "prob", probability = TRUE)
    predTest= attr(predTest, "probabilities")
    predTest2=prediction(predTest[,2], TEST$Outcome)
  }
  if (Var3==3){
    model <- glm(as.formula(paste('Outcome', paste(features, collapse=" + "), sep=" ~ ")), 
                 data = dt, 
                 family = binomial)
    predTrain <- predict(model, dt, type = "response")
    predTrain2=prediction(predTrain, dt$Outcome)
    predTest <- predict(model, TEST, type = "response")
    predTest2=prediction(predTest, TEST$Outcome)
  }
  
  if (Var3==4){
    model<- naiveBayes(Outcome ~ ., data = dt)
    predTrain <- predict(model, dt, type = "raw")
    predTrain2=prediction(predTrain[,2], dt$Outcome)
    predTest <- predict(model, TEST, type = "raw")
    predTest2=prediction(predTest[,2], TEST$Outcome)
  }
  
  AUC_TRAIN=performance(predTrain2,measure = "auc")@y.values[[1]]
  ROC_TRAIN=performance(predTrain2,"sens", "spec")
  YI=ROC_TRAIN@y.values[[1]]+ROC_TRAIN@x.values[[1]]-1
  ind_YI= match(max(YI),YI)
  THRESH= ROC_TRAIN@alpha.values[[1]][ind_YI]
  
  sensitivity_TRAIN=ROC_TRAIN@y.values[[1]][ind_YI]
  specificity_TRAIN=ROC_TRAIN@x.values[[1]][ind_YI]
  
  
  fscore_TRAIN_pr=performance(predTrain2,"prec", "rec")
  
  YI_pr=fscore_TRAIN_pr@y.values[[1]]+fscore_TRAIN_pr@x.values[[1]]-1
  ind_YI_pr= match(max(YI_pr[2:length(YI_pr)]),YI_pr)
  THRESH_pr=fscore_TRAIN_pr@alpha.values[[1]][ind_YI_pr]
  
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
  
  
  precision_TRAIN=(tp/(tp+fp))
  
  
    #prauc_TEST= pr.curve(TEST$Outcome,predTest)$auc.integral
    a=(performance(predTest2,"prec", "rec")@y.values[[1]])
    a=a[2:length(a)]
    b=diff(performance(predTest2,"prec", "rec")@x.values[[1]])
    prauc_TEST=sum(a*b)
  
    #prauc_TEST= pr.curve(TEST$Outcome,predTest[,2])$auc.integral
    
  
  AUC_TEST =performance( predTest2,measure = "auc")@y.values[[1]]
  ROC_TEST=performance( predTest2,"sens", "spec")
  YI=ROC_TEST@y.values[[1]]+ROC_TEST@x.values[[1]]-1
  ind_YI= match(max(YI),YI)
  THRESH=ROC_TEST@alpha.values[[1]][ind_YI]
  
  
  THRESH_TEST= ROC_TEST@alpha.values[[1]][ind_YI]
  
  
  fscore_TEST_pr=performance(predTest2,"prec", "rec")
  
  YI_pr=fscore_TEST_pr@y.values[[1]]+fscore_TEST_pr@x.values[[1]]-1
  ind_YI_pr= match(max(YI_pr[2:length(YI_pr)]),YI_pr)
  
  ind_YI_pr=length(which(as.numeric(predTest2@cutoffs[[1]])>0.5))
  
  tp=predTest2@tp[[1]][ind_YI_pr]
  tn=predTest2@tn[[1]][ind_YI_pr]
  fn=predTest2@fn[[1]][ind_YI_pr]
  fp=predTest2@fp[[1]][ind_YI_pr]
  
  
  
  
  if (tp==0){
    fscore_TEST=0
    recall=0
    precision=0
  }else{
    recall=tp/(tp+fn)
    precision=tp/(tp+fp)
  fscore_TEST=((1+(beta^2))*(precision*recall))/(((beta^2)*precision)+recall)
  
  }
  
  
  
  precision_TEST=precision
  recall_TEST=recall
  THRESH_TEST_pr= fscore_TEST_pr@alpha.values[[1]][ind_YI]
  
  if (Var3==3)
  {
    predTest_THRESH=predTest
  }else{
    predTest_THRESH=predTest[,2]
  }
  predTest_THRESH[which(predTest_THRESH>THRESH_TEST_pr)]=1 #THRESH_TEST
  predTest_THRESH[which(predTest_THRESH<=THRESH_TEST_pr)]=0
  Acc=Accuracy(predTest_THRESH, TEST$Outcome)
  
  sensitivity_TEST=ROC_TEST@y.values[[1]][ind_YI]
  specificity_TEST=ROC_TEST@x.values[[1]][ind_YI]
  
  return(c(AUC_TRAIN,AUC_TEST,sensitivity_TRAIN, sensitivity_TEST,specificity_TRAIN, specificity_TEST,THRESH,Acc,THRESH_TEST, fscore_TRAIN,fscore_TEST, precision_TEST, recall_TEST, prauc_TEST, tp,fp))
  
}