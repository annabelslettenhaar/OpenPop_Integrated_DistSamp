runMCMC_allcode <- function(model_setup, input_data, per_chain_info) {
  
  library(nimble)
  library(nimbleDistance)
  
  #assign("dHN", dHN, envir = .GlobalEnv)
  
  myModel <- nimbleModel(code = model_setup$modelCode,
                         data = input_data$nim.data, 
                         constants = input_data$nim.constants,
                         inits = per_chain_info$inits)

  CmyModel <- compileNimble(myModel)
  
  myMCMC <- buildMCMC(CmyModel, 
                      monitors = model_setup$modelParams)
  
  CmyMCMC <- compileNimble(myMCMC)
  
  results <- runMCMC(CmyMCMC, 
                     niter = model_setup$mcmcParams$niter, 
                     nburnin = model_setup$mcmcParams$nburn, 
                     thin = model_setup$mcmcParams$nthin, 
                     samplesAsCodaMCMC = TRUE, 
                     setSeed = per_chain_info$mySeed)
  
  return(results)
}
