Performance_Total <- function(Test_scanner,Classifier,Data,Data_scanner,Features) {
  
  R_data=Performance_TestData(Test_scanner,Data,Data_scanner)
  Set=R_data[[1]]
  Output_Set=R_data[[2]]
  TestSet=R_data[[3]]
  Output_TestSet=R_data[[4]]
  VTR=R_data[[5]]
  VTE=R_data[[6]]
  
  #Reset data with selected features
  Set=Set[Features]
  TestSet=TestSet[Features]
  
  ########
  
  Results_Classifier=sapply(Classifier,Performance_class,Set, Output_Set, TestSet, Output_TestSet,Features,VTR,VTE,Data)
  prauc_VAL=Results_Classifier[1,]
  prauc_TESTE=Results_Classifier[2,]
  
  return(list(prauc_VAL,prauc_TESTE,Features))
}