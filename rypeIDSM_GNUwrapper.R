library("tidyverse")
library("nimble")

rypeIDSM_GNUwrapper <- function(originSeed, chainSeed){
  
  ## Set seed for environment
  set.seed(originSeed)
  
  ## Set number of chains
  nchains <- 1
  
  ## Source relevant functions in "R" folder
  source("R/writeModelCode.R")
  source("R/setupModel.R")
  source("R/simulateInits.R")
  
  ## Load input data
  input_data <- readRDS("RypeData_forIM.rds")
  
  ## Load stored model switches/toggles
  toggles <- readRDS("ModelRunToggles.rds")
  
  ## Add definition-time if-else toggles to constants
  input_data$nim.constants$fitRodentCov <- toggles$fitRodentCov
  input_data$nim.constants$survVarT <- toggles$survVarT
  input_data$nim.constants$R_perF <- toggles$R_perF
  input_data$nim.constants$telemetryData <- toggles$telemetryData
  
  ## Write model code
  modelCode <- writeModelCode(survVarT = toggles$survVarT,
                              telemetryData = toggles$telemetryData)
  
  ## Setup for model using nimbleDistance::dHN
  model_setup <- setupModel(modelCode = modelCode,
                            R_perF = toggles$R_perF,
                            survVarT = toggles$survVarT, 
                            fitRodentCov = toggles$fitRodentCov,
                            nim.data = input_data$nim.data,
                            nim.constants = input_data$nim.constants,
                            testRun = toggles$testRun, 
                            nchains = nchains,
                            initVals.seed = chainSeed)
  ## Run model
  IDSM.out <- nimbleMCMC(code = model_setup$modelCode,
                         data = input_data$nim.data, 
                         constants = input_data$nim.constants,
                         inits = model_setup$initVals, 
                         monitors = model_setup$modelParams,
                         nchains = model_setup$mcmcParams$nchains, 
                         niter = model_setup$mcmcParams$niter, 
                         nburnin = model_setup$mcmcParams$nburn, 
                         thin = model_setup$mcmcParams$nthin, 
                         samplesAsCodaMCMC = TRUE, 
                         setSeed = chainSeed)

  ## Save output to RDS
  saveRDS(IDSM.out, 
          file = paste0("rypeIDSM_dHN_multiArea_realData_allAreas_GNU_", originSeed, "_", chainSeed, ".rds"))
  
}

args <- commandArgs(TRUE)

rypeIDSM_GNUwrapper(originSeed = as.integer(args[1]),
                    chainSeed = as.integer(args[2]))