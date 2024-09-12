storeToggles <- function(downloadData, R_perF, R_parent_drop0, sumR.Level,
                         survVarT, fitRodentCov, telemetryData, 
                         testRun, parallelMCMC){
  
  
  toggleList <- list(downloadData = downloadData,
                     R_perF = R_perF,
                     R_parent_drop0 = R_parent_drop0,
                     sumR.Level = sumR.Level,
                     survVarT = survVarT,
                     fitRodentCov = fitRodentCov,
                     telemetryData = telemetryData,
                     testRun = testRun,
                     parallelMCMC = parallelMCMC)
  
  saveRDS(toggleList, file = "ModelRunToggles.rds")
  
}