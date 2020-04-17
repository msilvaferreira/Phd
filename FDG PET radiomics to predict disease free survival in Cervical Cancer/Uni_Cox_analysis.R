Uni_Cox_analysis=function(Data, Bonferoni, treshold_p_value){
  
  # Description: This function performs a univariate cox regression analysis of the data frame Data. If demanded it also corrects 
  #for multiple testing using the Bonferoni method.
  #Input:
  #- Data: a data frame of the data to be analysed. Data should contain variables, survival time and survival outcome
  #- Bonferoni: If the user wants to correct for multiple testing Bonferoni should be the string 'yes'
  #treshold_p_value: p valye treshold to define statistical significance
  #Output:
  #- Returns a data frame with significant values (p value less or equal to treshold_p_value)
  
  DT=Data
  DT$Outcome <- NULL
  DT$Patient= NULL
  DT$Time= NULL
  
  
  
  survObj_Train <- Surv(time = Data$Time, event = Data$Outcome)
  
  ## Initialize univariate feature selection
  covariates <- colnames(DT)
  ##save each covariate as individial model formula
  univ_formulas <-
    sapply(covariates, function(x)
      as.formula(paste('survObj_Train ~', x))) 
  ##create a cox-model for each individual feature
  univ_models <-
    lapply(univ_formulas, function(x) {
      coxph(x, data = DT)
    })
  
  ## Calculate significance of each feature using the Wald test
  univ_results <- lapply(univ_models, function(x) {
    x <- summary(x)
    p.value <- signif(x$wald["pvalue"], digits = 2)
    beta <- signif(x$coef[1], digits = 2)
    #coeficient beta
    HR <- signif(x$coef[2], digits = 5)
    #Hasard ratios confidence Intervals
    HR.confint.lower <- signif(x$conf.int[, "lower .95"], 2)
    HR.confint.upper <- signif(x$conf.int[, "upper .95"], 2)
    wald.test <- paste0(signif(x$wald["test"], digits = 2),
                        " (",
                        HR.confint.lower,
                        "-",
                        HR.confint.upper,
                        ")")
    res <- c(beta, HR, wald.test, p.value) #save results
    names(res) <- c("beta", "HR", "wald.test",
                    "p.value")
    return(res)
  })
  
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  
  univariate_results <- as.data.frame(res)
  
  # Bonferoni correction
  if (Bonferoni=='yes'){
    
    alfa=treshold_p_value/(dim(univariate_results)[1])
    
    filtered_uv_results <-
      univariate_results[as.numeric(as.character(univariate_results$p.value)) <  alfa, ]
    
  }else{
    
    ##remove all features that do not have a (highly) significant according to p-value threshold
    filtered_uv_results <-
      univariate_results[as.numeric(as.character(univariate_results$p.value)) <  treshold_p_value, ]
    
    
  }
  return(filtered_uv_results)
}