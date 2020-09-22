HB_correction=function(Res_fc,alfa, Res1){
  
  # Description: This function performs a Holm Bonferoni correction
  
  #Output:
  #- Returns same data frame as initially with added column containing corrected p values (p_HB) 
  #and stop_c which returns the first feature index of the first non significant feature
  
  s=seq(1:(dim(Res_fc)[1]))
  
  Res1=Res_fc[order(Res_fc$p.value),]
  Res1['seq']=s
  Res1['p_HB']=alfa/((dim(Res1)[1])-Res1$seq+1)
  
  
  Stop_c=which(as.numeric(as.vector(Res1$p.value))>as.numeric(as.vector(Res1$p_HB)))[1]
  
  return(Stop_c, Res1)

}