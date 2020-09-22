Transform_ComBat= function (DATA_radiomics_all, Data_Patient, Data_scanner){
  
  source("C:\\Users\\f\\Documents\\Phd\\Project 1\\code\\Paper 1\\combat.R")
  #This function applies necessary transformation to data frame to apply combate
  #Input:
  #DATA_radiomics_all= Data frame with radiomics features
  #Data_Patient = data frame with patient list
  #Data_scanner = data frame with information regarding scanner for each patient
  #Output:
  #DATA_radiomics_all = data frame with features harmonised with comBat

DATA_radiomics_all$Patient=NULL
Data_Combat <- transpose(DATA_radiomics_all)

colnames(Data_Combat) <- rownames(DATA_radiomics_all)
rownames(Data_Combat) <- colnames(DATA_radiomics_all)


# Transform scanners into numbers

Scanner=Data_scanner$Scanner

Scanner_cat <- matrix(0, length(Scanner), 1)

for (j in 1:(length(Scanner)))
{
  
  
  if (Scanner[j]=='A'){
       Scanner_cat[j]=1
  }
  if (Scanner[j]=='B'){
    Scanner_cat[j]=2
  }
  
  if (Scanner[j]=='EF'){
    Scanner_cat[j]=3
  }
  if (Scanner[j]=='G'){
    Scanner_cat[j]=4
  }
  if (Scanner[j]=='H'){
    Scanner_cat[j]=5
  }
}

batch =c(Scanner_cat)


data_harmonized <- combat(dat=Data_Combat, batch=batch)

dt_harm=data_harmonized$dat.combat

Data_again<- as.data.frame(t(dt_harm))

colnames(Data_again) <- colnames(DATA_radiomics_all)
rownames(Data_again) <- rownames(DATA_radiomics_all)

DATA_radiomics_all=Data_again
DATA_radiomics_all['Patient']=Data_Patient$Patient

return(DATA_radiomics_all)
}