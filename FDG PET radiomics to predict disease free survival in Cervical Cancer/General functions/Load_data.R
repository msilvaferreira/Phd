Load_data= function (General_path,data_type, Data_Patient, Disccretisation_method, interpolation){
  
  
  #  Load_data is a function which returns the OR radiomics, TLR radiomics, clinical and treatment data, or Output(DFS)
  # Input:
  # data_type= string with data_type ("clinical", "Output","treatment", "OR_Radiomics" or "TLR_Radiomics"
  # Data_Patient = Data frame with list of patients
  # Disccretisation_method 'FBW' or 'FBN'. 'None' if any of these are applicable
  # Output:
  # Data= data frame with patient list(first column) and requested data type
  # 
  
  
  
  
  #Clinical Data
  if (data_type=="clinical"){
    DATA_C <- read.csv(paste(General_path,'Excell_data\\Other\\Data_Clinical_169patient.csv',sep=""), sep = ';')
    Data=DATA_C
  }
  if (data_type=="treatment"){
    DATA_C <- read.csv(paste(General_path,'Excell_data\\Other\\patient_treatment.csv',sep=""), sep = ';')
    Data=DATA_C
  }
  if (data_type=="Output"){
    DATA_output <- read.csv(paste(General_path,'Excell_data\\Other\\Data_DFS_169patient.csv',sep=""), sep = ',')
    DATA_output$Institute=NULL
    Data=DATA_output
  }
  if (data_type=="Scanner"){
    DATA_output <- read.csv(paste(General_path,'Excell_data\\Other\\Data169_scanner_binary.csv',sep=""), sep = ',')
    Data=DATA_output
  }
  
  
  if (data_type=="OR_Radiomics" || data_type=="TLR_Radiomics"){
    if (Disccretisation_method=='FBW'){
      
    # Radiomics Features
    DATA_radiomics_Liege_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Liege\\data_radiomics_CHU_Liege_all.csv',sep=""), sep=",")
    DATA_radiomics_Barcelona_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Barcelona\\data_radiomics_Barcelona_all.csv',sep=""), sep=",")
    if (interpolation=='True'){
      DATA_radiomics_Brest_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Brest\\data_radiomics_Brest_FWB_all_interpolation.csv',sep=""), sep=",")
      DATA_radiomics_Mcguill_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\McGuill\\data_radiomics_McGuill_FWB_all_interpolation.csv',sep=""), sep=",")
      DATA_radiomics_Nantes_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Nantes\\Radiomics_PET_FBW_interpolation_linear_Nantes.csv',sep=""), sep=",")
      
    }else{
    DATA_radiomics_Brest_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Brest\\data_radiomics_Brest_all.csv',sep=""), sep=",")
    DATA_radiomics_Mcguill_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\McGuill\\data_radiomics_McGuill_all.csv',sep=""), sep=",")
    DATA_radiomics_Nantes_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Nantes\\Radiomics_PET_FBW_Nantes.csv',sep=""), sep=",")
    }
    }
    
    
    
    if (Disccretisation_method=='FBN'){
      
      # Radiomics Features
      DATA_radiomics_Liege_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Liege\\data_radiomics_Liege_FNB_all.csv',sep=""), sep=",")
      DATA_radiomics_Barcelona_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Barcelona\\data_radiomics_Barcelona_FNB_all.csv',sep=""), sep=",")
      if (interpolation=='True'){
        DATA_radiomics_Brest_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Brest\\Radiomics_PET_FBN_interpolation.csv',sep=""), sep=",")
        DATA_radiomics_Mcguill_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\McGuill\\Radiomics_PET_FBN_interpolation.csv',sep=""), sep=",")
        DATA_radiomics_Nantes_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Nantes\\Radiomics_PET_FBN_interpolation_linear_Nantes.csv',sep=""), sep=",")
        
      }else{
      DATA_radiomics_Brest_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Brest\\data_radiomics_Brest_FNB_all.csv',sep=""), sep=",")
      DATA_radiomics_Mcguill_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\McGuill\\data_radiomics_McGuill_FNB_all.csv',sep=""), sep=",")
      DATA_radiomics_Nantes_FBW <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Nantes\\Radiomics_PET_FBN_Nantes.csv',sep=""), sep=",")
      }
    }
      # Bind radiomics from the different institutes
      DATA_radiomics_all <- rbind(DATA_radiomics_Liege_FBW,  DATA_radiomics_Barcelona_FBW )
      DATA_radiomics_all <- rbind(DATA_radiomics_all, DATA_radiomics_Mcguill_FBW)
      DATA_radiomics_all <- rbind(DATA_radiomics_all, DATA_radiomics_Brest_FBW )
      DATA_radiomics_all <- rbind(DATA_radiomics_all, DATA_radiomics_Nantes_FBW )
    
    
    DATA_radiomics_all=merge( DATA_radiomics_all,Data_Patient,by="Patient")
    
    Data=DATA_radiomics_all
    
    
  }
  
  
  
  if (data_type=="TLR_Radiomics"){
    
    if (Disccretisation_method=='FBW'){
    
    # Radiomics Features
    DATA_radiomics_Liege_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Liege\\data_radiomics_CHU_Liege_FWB_all_liver.csv',sep=""), sep=",")
    DATA_radiomics_Barcelona_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Barcelona\\data_radiomics_Barcelona_FWB_all_liver.csv',sep=""), sep=",")
    if (interpolation=='True'){
      DATA_radiomics_Brest_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Brest\\data_radiomics_CHRU_FWB_all_liver_interpolation.csv',sep=""), sep=",")
      DATA_radiomics_Mcguill_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\McGuill\\data_radiomics_McGuill_FWB_all_liver_interpolation.csv',sep=""), sep=",")
      DATA_radiomics_Nantes_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Nantes\\Radiomics_PET_FBW_interpolation_linear_liver_Nantes.csv',sep=""), sep=",")
      
    }else{
    DATA_radiomics_Brest_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Brest\\data_radiomics_CHRU_FWB_all_liver.csv',sep=""), sep=",")
    DATA_radiomics_Mcguill_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\McGuill\\data_radiomics_McGuill_FWB_all_liver.csv',sep=""), sep=",")
    DATA_radiomics_Nantes_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Nantes\\Radiomics_PET_FBW_liver_Nantes.csv',sep=""), sep=",")
    }
    }
    
    if (Disccretisation_method=='FBN'){
      DATA_radiomics_Liege_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Liege\\data_radiomics_CHU_Liege_FNB_all_liver.csv',sep=""), sep=",")
      DATA_radiomics_Barcelona_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Barcelona\\data_radiomics_Barcelona_FNB_all_liver.csv',sep=""), sep=",")
      if (interpolation=='True'){
        DATA_radiomics_Brest_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Brest\\Radiomics_PET_FBN_interpolation_liver.csv',sep=""), sep=",")
        DATA_radiomics_Mcguill_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\McGuill\\Radiomics_PET_FBN_interpolation_liver.csv',sep=""), sep=",")
        DATA_radiomics_Nantes_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Nantes\\Radiomics_PET_FBN_interpolation_linear_liver_Nantes.csv',sep=""), sep=",")
        
      }else{
      DATA_radiomics_Brest_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Brest\\data_radiomics_CHRU_FNB_all_liver.csv',sep=""), sep=",")
      DATA_radiomics_Mcguill_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\McGuill\\data_radiomics_McGuill_FNB_all_liver.csv',sep=""), sep=",")
      DATA_radiomics_Nantes_FBW_liver <- read.csv(paste(General_path,'Excell_data\\Radiomics features\\Nantes\\Radiomics_PET_FBN_liver_Nantes.csv',sep=""), sep=",")
      }
    }
    
    # Bind radiomics from the different institutes
    DATA_radiomics_all_liver <- rbind(DATA_radiomics_Liege_FBW_liver,  DATA_radiomics_Barcelona_FBW_liver )
    DATA_radiomics_all_liver <- rbind(DATA_radiomics_all_liver, DATA_radiomics_Mcguill_FBW_liver)
    DATA_radiomics_all_liver <- rbind(DATA_radiomics_all_liver, DATA_radiomics_Brest_FBW_liver )
    DATA_radiomics_all_liver <- rbind(DATA_radiomics_all_liver, DATA_radiomics_Nantes_FBW_liver )
    
    DATA_radiomics_all_liver=merge( DATA_radiomics_all_liver,Data_Patient,by="Patient")
    
    old_DATA_radiomics_all_liver=DATA_radiomics_all_liver
    
    #Ratio features
    DATA_radiomics_all_liver[2:dim(DATA_radiomics_all)[2]]=DATA_radiomics_all[2:dim(DATA_radiomics_all)[2]]/DATA_radiomics_all_liver[2:dim(DATA_radiomics_all)[2]]
    
    
    #Replace by OR features features related with shape
    for (i in 173:196){
      DATA_radiomics_all_liver[i]=old_DATA_radiomics_all_liver[i]
      
    }
    
    Data=DATA_radiomics_all_liver
  }
  return(Data)
}