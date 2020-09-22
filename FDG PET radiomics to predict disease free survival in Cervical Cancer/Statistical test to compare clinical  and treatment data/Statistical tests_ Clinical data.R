################
library(data.table)
source("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Code\\Paper 1\\Load_data.R")

#Construct tables for variables outcome for 4 groups
Risk_matrix <- function(car,car_list,Scanner_list,Data){
  Performance <-matrix(c(1:(length(Scanner_list)))*0, nrow = length(car_list),ncol = length(Scanner_list))
  
  
  for (i in 1:(dim(Data))[1])
  {
    r=which(car_list==Data[car][[1]][i])
    c=which(Scanner_list==Data['Scanner'][[1]][i])
    
    
    Performance[r,c]<- Performance[r,c]+1
  }
      

  return(Performance)
  
}


#Load data

Data_Patient= read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Excell data\\Data_169patient.csv", sep = ';')
Data_Patient= data.frame("Patient" = Data_Patient$Patient)

DATA_C=Load_data("clinical", Data_Patient, 'None','None')
DATA_C$Institute=NULL
DATA_C=merge(DATA_C,Data_Patient, by='Patient')

DATA_treatment=Load_data("treatment", Data_Patient, 'None','None')
DATA_treatment=merge(DATA_treatment,Data_Patient, by='Patient')

Data_scanner=read.csv("C:\\Users\\Marta\\Google Drive\\OneDrive - student.uliege.be\\PHD\\Data\\Excell data\\Data169_scanner.csv", sep = ';')
Data_scanner=merge(Data_scanner,Data_Patient,by='Patient')

D_enter=merge(DATA_C,Data_scanner,by='Patient')
D_enter=merge(D_enter,DATA_treatment,by='Patient')


####################
#FIGO
Performance_FIGO = Risk_matrix("FIGO",c(1,2,3,4),c("A","EF","H"),D_enter)
TABLE_FIGO=as.table(Performance_FIGO)
rownames(TABLE_FIGO) <- c("1", "2",'3','4')
colnames(TABLE_FIGO) <- c("A","EF","H")

FIGO_STATS= chisq.test(Performance_FIGO, correct =FALSE)

###################################

# Histology
Performance_Histology = Risk_matrix("Histology",c(1,0),c("A","EF","H"),D_enter)
TABLE_Histology=as.table(Performance_Histology)
rownames(TABLE_Histology) <- c("SSC", "Other")
colnames(TABLE_Histology) <- c("A","EF","H")

Histology_STATS= chisq.test(Performance_Histology, correct =FALSE)

###################################

#Metastasis
Performance_Metastasis = Risk_matrix("Metastasis",c(1,0),c("A","EF","H"),D_enter)
TABLE_Metastasis=as.table(Performance_Metastasis)
rownames(TABLE_Metastasis) <- c("Metastasis", "No Metastasis")
colnames(TABLE_Metastasis) <- c("A","EF","H")

Metastasis_STATS= chisq.test(Performance_Metastasis, correct =FALSE)

###########################
#Age

dANOVA=data.frame('Age'=D_enter$Age, 'Scanner'=D_enter$Scanner)
table_Age_df=data.frame('A'=mean(D_enter[which(D_enter=="A"),]$Age),'EF'=mean(D_enter[which(D_enter=="EF"),]$Age),'H'=mean(D_enter[which(D_enter=="H"),]$Age))
table_Age=as.data.table(table_Age_df)
rownames(table_Age)=c('Age')

AGE_STATS= aov(Age ~ Scanner, data=dANOVA)
summary(AGE_STATS)

###################################

#type_RTE
Performance_type_RTE = Risk_matrix("type_RTE",c(1,0),c("A","EF","H"),D_enter)
TABLE_type_RTE=as.table(Performance_type_RTE)
rownames(TABLE_type_RTE) <- c("3D", "No 3D")
colnames(TABLE_type_RTE) <- c("A","EF","H")

type_RTE_STATS= chisq.test(Performance_type_RTE, correct =FALSE)

###################################

#dose_ggl
Performance_dose_ggl = Risk_matrix("dose_ggl",c(1,0),c("A","EF","H"),D_enter)
TABLE_dose_ggl=as.table(Performance_dose_ggl)
rownames(TABLE_dose_ggl) <- c("Dose in ggl", "No dose in ggl")
colnames(TABLE_dose_ggl) <- c("A","EF","H")

dose_ggl_STATS= chisq.test(Performance_dose_ggl, correct =FALSE)

####################
#type_BrachRad
Performance_type_BrachRad = Risk_matrix("type_BrachRad",c(1,2,3,4),c("A","EF","H"),D_enter)
TABLE_type_BrachRad=as.table(Performance_type_BrachRad)
rownames(TABLE_type_BrachRad) <- c("1", "2",'3','4')
colnames(TABLE_type_BrachRad) <- c("A","EF","H")

type_BrachRad_STATS= chisq.test(Performance_type_BrachRad, correct =FALSE)
