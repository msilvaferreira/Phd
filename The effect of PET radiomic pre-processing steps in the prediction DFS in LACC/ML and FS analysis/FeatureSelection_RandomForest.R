###########################################
## Feature selection with Random Forest ###
###########################################

# Function selects features through the embed randomForest feature selection
#Input:
#TR: features data frame
#Z= seed value
#Outcome= vector of Output values
#Criteria= string with criteria of feature importance ('MeanDecreaseAccuracy' or 'MeanDecreaseGini')
#NR_Features= Number of features
#Output:
#TR_Return= data frame with reduced number of features

FeatureSelection_RandomForest <- function(TR,Z,Outcome,criteria,NR_Features){
  TR['Outcome']=Outcome
  set.seed(Z)
  m <- randomForest(Outcome ~ ., data = TR, importance = TRUE)
  imp= as.data.frame(importance(m)) #change according to model (model1 or model2)
  newdata_imp <- imp[order(-imp[criteria]),]
  hd_newdata_imp=row.names(newdata_imp)
  new_hd_newdata_imp=hd_newdata_imp[1:NR_Features]
  TR['Outcome']=NULL
  TR_Return=TR[c(new_hd_newdata_imp)]
  
  return(TR_Return)
  
}

