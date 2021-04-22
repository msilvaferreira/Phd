FeatureSelection_MLMR <- function(Data,Output,NR_Features,direction,method){
  
  #Score function
  evaluate_score <- function(remaining,selection, Data,Output,method){
    if (length(selection)==0){
      score=abs(cor(Data[remaining],Output, use='all.obs', method=method ))
      score[is.na(score)] <- 0
    }else{
      s1=abs(cor(Data[remaining], Output, use='all.obs', method=method ))
      s1[is.na(s1)] <- 0
      s2=(sum(abs(cor(Data[remaining],Data[selection], use='all.obs', method=method )))/dim(Data[selection])[2])
      s2[is.na(s2)] <- 0
      score=s1-s2
      
    }
    
    return(score[1,1])
  }
  
  
  Output=ifelse((as.numeric(Output)) > 1, 1, 0)
  
  #Forward selection
  if (direction=='up'){
    
    forward = function(remaining_initial,selection_initial,Data, steps, Output,method){
      remaining=remaining_initial
      selection=selection_initial
      while(length(selection) < steps){
        
        scores=sapply(remaining,evaluate_score, selection, Data,Output,method)
        indexsel = which.max(scores)
        indexsel=indexsel[1]
        selection =c(remaining[indexsel], selection)
        remaining = remaining[-indexsel] 
      }
      
      return(c(selection,scores[indexsel]))
      
    }  
    
    remaining_initial = colnames(Data) #Unselected players
    selection_initial = c()
    
    Results=forward(remaining_initial,selection_initial,Data,NR_Features, Output,method)
  }
  
  #Backward Selection
  if (direction=='down'){
    
    backward = function(remaining_initial,selection_initial,Data, steps, Output,method){
      remaining=remaining_initial
      selection=selection_initial
      while(length(selection) > steps){
        
        
        scores=sapply(selection,evaluate_score, remaining, Data,Output,method)
        indexsel = which.min(scores)
        indexsel=indexsel[1]
        remaining =c(selection[indexsel], remaining)
        selection = selection[-indexsel]
      }
      
      return(c(selection,scores[indexsel]))
      
    } 
    
    
    remaining_initial = c() 
    selection_initial = colnames(Data)
    
    Results=backward(remaining_initial,selection_initial,Data,NR_Features, Output,method)
  }
  
  
  return(Data[Results[1:(length(Results)-1)]])
  
  
}