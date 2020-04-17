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


my_RF_to_BOOTSTRAP_Sensitivity_AUC <- function(data, indices, TEST,Var3,features){
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
  
  AUC_TEST =performance( predTest2,measure = "auc")@y.values[[1]]
  ROC_TEST=performance( predTest2,"sens", "spec")
  YI=ROC_TEST@y.values[[1]]+ROC_TEST@x.values[[1]]-1
  ind_YI= match(max(YI),YI)
  
  THRESH_TEST= ROC_TEST@alpha.values[[1]][ind_YI]
  if (Var3==3)
  {
    predTest_THRESH=predTest
  }else{
    predTest_THRESH=predTest[,2]
  }
  predTest_THRESH[which(predTest_THRESH>THRESH_TEST)]=1
  predTest_THRESH[which(predTest_THRESH<=THRESH_TEST)]=0
  Acc=Accuracy(predTest_THRESH, TEST$Outcome)
  
  sensitivity_TEST=ROC_TEST@y.values[[1]][ind_YI]
  specificity_TEST=ROC_TEST@x.values[[1]][ind_YI]
  
  return(c(AUC_TRAIN,AUC_TEST,sensitivity_TRAIN, sensitivity_TEST,specificity_TRAIN, specificity_TEST,THRESH,Acc,THRESH_TEST))
  
}