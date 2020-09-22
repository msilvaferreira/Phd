rm(list = ls())

##### libraries ######

library(survival)
library(survminer)
library(standardize)
library(caret)
library(dplyr)
library(pROC)
##### load functions #####

General_path='C:\\Users\\f\\Documents\\Phd\\Project 1\\'

source(paste(General_path,"code\\Paper 1\\Uni_Cox_analysis.R",sep=""))
source(paste(General_path,"code\\Paper 1\\Load_data.R",sep=""))
source(paste(General_path,"code\\Paper 1\\Transform_ComBat.R",sep=""))
source(paste(General_path,"code\\Paper 1\\HB_correction.R",sep=""))
#### read data ######

Data_Patient= read.csv(paste(General_path,"Excell_data\\Other\\Data_140patient.csv",sep=""), sep = ';')
Data_Patient= data.frame("Patient" = Data_Patient$Patient)

DATA_C=Load_data(General_path,"clinical", Data_Patient, 'None','None')
DATA_C$Institute=NULL

DATA_treatment=Load_data(General_path,"treatment", Data_Patient, 'None','None')
DATA_treatment=merge(DATA_treatment,Data_Patient, by='Patient')

Data_scanner=Load_data(General_path,"Scanner", Data_Patient, 'None','None') 
Data_scanner=merge(Data_scanner,Data_Patient,by='Patient')

DATA_output=Load_data(General_path,"Output", Data_Patient, 'None','None')

Data_Radiomics=Load_data(General_path,"OR_Radiomics", Data_Patient, 'FBN','False')

Data_TLR_Radiomics=Load_data(General_path,"TLR_Radiomics",Data_Patient, 'FBW','True')

Data_Liver_Radiomics=Load_data(General_path,"Liver_Radiomics",Data_Patient, 'FBN','False')


### merge features with clinics ###

Data =merge( DATA_C,DATA_output,by="Patient")
Data =merge( Data,DATA_treatment,by="Patient")
Data =merge( Data,Data_scanner,by="Patient")
Data =merge( Data,Data_Radiomics,by="Patient")
Data =merge( Data,Data_TLR_Radiomics,by="Patient")
Data =merge( Data,Data_Liver_Radiomics,by="Patient")

D_T=Data
####################

D_T['Outcome']=Data$DFS
D_T['Patient']=Data$Patient
D_T['Outcome']=as.factor(D_T$Outcome)


#Split into train and test according to test scanner
Z=1
ratio=0.8
set.seed(Z)

D_A=Data_scanner[which(Data_scanner$Scanner=='A'),]

D_EF=Data_scanner[which(Data_scanner$Scanner=='EF'),]

set.seed(Z)
train_A <- sample(nrow(D_A), ratio*nrow(D_A), replace = FALSE)#
Set_A <- D_A[train_A,]
TestSet_A <- D_A[-train_A,]

set.seed(Z)
train_EF <- sample(nrow(D_EF), ratio*nrow(D_EF), replace = FALSE)#
Set_EF <- D_EF[train_EF,]
TestSet_EF <- D_EF[-train_EF,]

Set=rbind(Set_A,Set_EF)
Set$Institute=NULL
Set$Scanner=NULL

Set=merge(Set,D_T,by='Patient')
Set$Patient <- NULL


TestSet=rbind(TestSet_A,TestSet_EF)
TestSet$Institute=NULL
TestSet$Scanner=NULL

TestSet=merge(TestSet,D_T,by='Patient')
TestSet$Patient <- NULL

#################################
Data=Set
Data$Outcome=NULL


Data$Scanner=NULL
Data$Institute=NULL

Data_original=Data

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

#Remove correlated features
New_Res1=Res
df1=Data[,rownames(New_Res1)]
df2 = cor(df1)
hc = findCorrelation(df2, cutoff=0.9) # putt any value as a "cutoff" 
hc = sort(hc)
reduced_Data = df1[,-c(hc)]
Data_uncorr=data.frame(reduced_Data)

Res_fc=Res
Res_fc=Res_fc[colnames(Data_uncorr),]

#Correct p-values with HB method

alfa=0.05
Res_HB=HB_correction(Res_fc,alfa, Res1)
Stop_c=Res_HB[1]
Res1=Res_HB[2]
  
f_after=c( 'Histology', 'IH_energy_64')

####

#Multivariate statistics

#Data_uncorr=Data

survObj_Train <- Surv(time = Time_Data, event = Outcome_Data)
measurevar <- "survObj_Train"
#groupvars  <- colnames(reduced_Data)
#groupvars  <- colnames(Data_uncorr)
groupvars= f_after

#groupvars=c('IH_maxGrad_64','IVH_RVRI_40_32')

for (i in 1:length( groupvars))
{
  if (i==1){
    
    st=toString(groupvars[i])
    
  }
  if (i>1){
    st= paste(c(st,toString(groupvars[i])),collapse= '+')
    
  }
  
}
coxRegMulti <- coxph(as.formula(paste(measurevar, st, sep=" ~ ")), data = Data_uncorr, x = TRUE, y = TRUE)
r1=summary(coxRegMulti)
r2=r1$conf.int
r3=r1$coefficients

beta=matrix(0, (dim(r2)[1]), 1)
hazard=matrix(0, (dim(r2)[1]), 1)
hazardlower=matrix(0, (dim(r2)[1]), 1)
hazardup=matrix(0, (dim(r2)[1]), 1)
pvalue=matrix(0, (dim(r2)[1]), 1)

for (i in 1:(dim(r2)[1]))
{
  beta[i]=r3[,1][[i]]
  hazard[i]=r3[,2][[i]]
  hazardlower[i]=r2[,3][[i]]
  hazardup[i]=r2[,4][[i]]
  pvalue[i]=r3[,5][[i]]
  
}
#write.csv(multi_dataframe, file = "Multivariate_Analysis_TLR_FBW.csv")
multi_dataframe=data.frame('Feature'= groupvars,'Beta'=beta, 'Hazard'=hazard, 'Hazard lower'=hazardlower, 'Hazard Up'=hazardup, 'P value'=pvalue)

#Other p value correction method

s=seq(1:(dim(multi_dataframe)[1]))

Res1=multi_dataframe[order(multi_dataframe$P.value),]
Res1['seq']=s
alfa=0.05
Res1['p_HB']=alfa/(dim(Res1)[1]-Res1$seq+1)
Res1['p_B']=alfa/(dim(Res1)[1])

Stop_c=which(as.numeric(as.vector(Res1$P.value))>as.numeric(as.vector(Res1$p_HB)))[1]

#write.csv(Res1, file = "Muti_Analysis.csv")

feature=c('GLDZM_LDE_0.1')
DT=Data_original[feature]

roc = roc(Outcome_Data,DT$GLDZM_LDE_0.1)

plot(roc)
YI_m=abs(roc$sensitivities+roc$specificities-1)
ind_YI_m= match(max(YI_m),YI_m)
THRESH_m= roc$thresholds[ind_YI_m]

roc$direction

### Test ####

DT=TestSet[feature]

Outcome_Data=TestSet$DFS
Time_Data=TestSet$Time


roc = roc(Outcome_Data,DT$GLDZM_LDE_0.1, direction="<")

plot(roc)

#######

#KM analysis

Data_KM= DT
#Data_KM= TestSet


for (i in 1:(dim(Data_KM))[1])
{
  
  if ((Data_KM$GLDZM_LDE_0.1)[i]>THRESH_m){
    
    Data_KM$GLDZM_LDE_0.1[i]=1
  }
  else {
    Data_KM$GLDZM_LDE_0.1[i]=0
  }
}

survObj<- Surv(time = Time_Data, event = Outcome_Data)

# kaplan meyer curve (change km variable)
library(survminer)
fit1 <- survfit(survObj~ Data_KM$GLDZM_LDE_0.1, data = Data_KM)
ggsurvplot(fit1, data =Data_KM, pval = TRUE)

#Analise Tp, fp, etc

tp=length(which((Data_KM$GLDZM_LDE_0.1==1  & Outcome_Data==1)==TRUE))
tn=length(which((Data_KM$GLDZM_LDE_0.1==0  & Outcome_Data==0)==TRUE))
fp=length(which((Data_KM$GLDZM_LDE_0.1==1  & Outcome_Data==0)==TRUE))
fn=length(which((Data_KM$GLDZM_LDE_0.1==0  & Outcome_Data==1)==TRUE))

recall=tp/(tp+fn)
precision=tp/(tp+fp)
beta=1
fscore=((1+(beta^2))*(precision*recall))/(((beta^2)*precision)+recall)

#write.csv(Res1, file = "Multi.csv")
#215+((215-90)*3)+((215-90)*2)
#215+((215-90)*3)+((215-90)*2)- (24*6)
#8