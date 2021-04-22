chose_descretisation_nr <- function(Disccretisation,DATA_radiomics_all) {
  
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
  
  if (Disccretisation=='All'){
    
    DATA_radiomics_all= DATA_radiomics_all
    
    
  }
  
  return(DATA_radiomics_all)
}