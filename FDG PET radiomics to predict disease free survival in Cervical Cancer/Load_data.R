Load_data= function (data_type){
  
  # Load_data is a function which returns the OR radiomics, TLR radiomics, clinical data or Output(DFS)
  #Input:
  #data_type= string with data_type ("clinical", "Output", "OR_Radiomics" or "TLR_Radiomics"
  #Output:
  #Data= data frame with patient list(first column) and requested data type
  
  Data_Patient= read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Excell data\\Data_148patient.csv", sep = ';')
  Data_Patient= data.frame("Patient" = Data_Patient$Patient)
  
  #Clinical Data
  if (data_type=="clinical"){
    DATA_C <- read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Data_Clinical_Discrete_secondHist.csv", sep = ';')
    Data=DATA_C
  }
  if (data_type=="Output"){
    DATA_output <- read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Data_DFS.csv", sep = ';')
    DATA_output$Institute=NULL
    Data=DATA_output
  }
  
  
  if (data_type=="OR_Radiomics" || data_type=="TLR_Radiomics"){
    # Radiomics Features
    DATA_radiomics_Liege_FBW <- read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Radiomics features_Liege_Barcelona\\CHU_Liege\\data_radiomics_CHU_Liege_all.csv", sep=";")
    DATA_radiomics_Barcelona_FBW <- read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Radiomics features_Liege_Barcelona\\Barcelona\\data_radiomics_Barcelona_all.csv", sep=";")
    DATA_radiomics_Brest_FBW <- read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Radiomics features_Liege_Barcelona\\CHRU\\data_radiomics_Brest_all.csv", sep=";")
    DATA_radiomics_Mcguill_FBW <- read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Radiomics features_Liege_Barcelona\\McGuill\\data_radiomics_McGuill_all.csv", sep=";")
    
    # Bind radiomics from the different institutes
    DATA_radiomics_all <- rbind(DATA_radiomics_Liege_FBW,  DATA_radiomics_Barcelona_FBW )
    DATA_radiomics_all <- rbind(DATA_radiomics_all, DATA_radiomics_Mcguill_FBW)
    DATA_radiomics_all <- rbind(DATA_radiomics_all, DATA_radiomics_Brest_FBW )
    
    DATA_radiomics_all=merge( DATA_radiomics_all,Data_Patient,by="Patient")
    
    Data=DATA_radiomics_all
    
  }
  
  if (data_type=="TLR_Radiomics"){
    
    # Radiomics Features
    DATA_radiomics_Liege_FBW_liver <- read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Radiomics features_Liege_Barcelona\\CHU_Liege\\data_radiomics_CHU_Liege_FWB_all_liver.csv", sep=";")
    DATA_radiomics_Barcelona_FBW_liver <- read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Radiomics features_Liege_Barcelona\\Barcelona\\data_radiomics_Barcelona_FWB_all_liver.csv", sep=";")
    DATA_radiomics_Brest_FBW_liver <- read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Radiomics features_Liege_Barcelona\\CHRU\\data_radiomics_CHRU_fixbinwidth_all_interpolation_nearest_liver.csv", sep=";")
    DATA_radiomics_Mcguill_FBW_liver <- read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Radiomics features_Liege_Barcelona\\McGuill\\data_radiomics_McGuill_fixbinwidth_all_interpolation_nearest_liver.csv", sep=";")
    
    # Bind radiomics from the different institutes
    DATA_radiomics_all_liver <- rbind(DATA_radiomics_Liege_FBW_liver,  DATA_radiomics_Barcelona_FBW_liver )
    DATA_radiomics_all_liver <- rbind(DATA_radiomics_all_liver, DATA_radiomics_Mcguill_FBW_liver)
    DATA_radiomics_all_liver <- rbind(DATA_radiomics_all_liver, DATA_radiomics_Brest_FBW_liver )
    
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