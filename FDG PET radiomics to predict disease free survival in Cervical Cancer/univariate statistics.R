rm(list = ls())

##### libraries ######

library(survival)
library(survminer)
library(standardize)

##### load functions #####

source("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Code\\Paper 1\\Uni_Cox_analysis.R")
source("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Code\\Paper 1\\Load_data.R")

#### read data ######

Data_Patient= read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Excell data\\Data_148patient.csv", sep = ';')
Data_Patient= data.frame("Patient" = Data_Patient$Patient)

DATA_C=Load_data("clinical")
DATA_C$Institute=NULL

DATA_output=Load_data("Output")

Data_Radiomics=Load_data("OR_Radiomics")

Data_TLR_Radiomics=Load_data("TLR_Radiomics")

### merge features with clinics ###

Data =merge( DATA_C,DATA_output,by="Patient")
Data =merge( Data,Data_Radiomics,by="Patient")
Data =merge( Data,Data_TLR_Radiomics,by="Patient")

Outcome_Data=Data$DFS
Data$DFS=NULL
Data$Patient=NULL
Time_Data=Data$Time
Data$Time=NULL

#Standardize data
Data=data.frame(scale(Data))

Data=Data[is.finite(colSums(Data)) ]
Data=Data[, sapply(Data, var) != 0]

Data['Time']=Time_Data
Data['Outcome']=Outcome_Data

Res= Uni_Cox_analysis(Data,'no',0.05)


