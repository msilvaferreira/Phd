Discretise=function(DATA_radiomics_all,Disccretisation){
  
  # Description: This function selects features according to desiderd descretisation
  #Input:
  #DATA_radiomics_all- data frame with features
  #Disccretisation-desiderd descretisation
  #Output:
  #- Returns same data frame as initially with added column containing corrected p values (p_HB) 
  #and stop_c which returns the first feature index of the first non significant feature

if (Disccretisation==0.5){
  #0.5
  DATA_radiomics_all= DATA_radiomics_all[-c(2:646)] # discretisation 0.5
  
  
}

if (Disccretisation==0.2){
  #0.2
  DATA_radiomics_all= DATA_radiomics_all[-c(2:431)] # discretisation 0.2
  DATA_radiomics_all= DATA_radiomics_all[-c(217:431)] # discretisation 0.2
  
  
}

if (Disccretisation==0.1){
  #0.1
  DATA_radiomics_all= DATA_radiomics_all[c(1:431)] # discretisation 0.1
  DATA_radiomics_all= DATA_radiomics_all[-c(2:216)] # discretisation 0.1
  
  
  
}

if (Disccretisation==0.05){
  #0.05
  DATA_radiomics_all= DATA_radiomics_all[-c(217:861)] # discretisation 0.05
  
  
}

if (Disccretisation==32){
  
  DATA_radiomics_all= DATA_radiomics_all[-c(217:431)]
  
  
}

if (Disccretisation==64){
  
  DATA_radiomics_all= DATA_radiomics_all[-c(2:216)]
  
  
}
  
return(DATA_radiomics_all)
}