#' Rewrite results served as RDS as Excel/CSV tables
#' 
#' Reads the RDS files containing posterior summaries of parameter estimates, 
#' calculated generation times, and variance decompositions and saves them as
#' .xlsx / .csv for publishing. 
#'
#' @return
#' @export
#'
#' @examples

writeResultsTables <- function(){
  postSum_main <- readRDS("PosteriorSummaries_byArea.rds")
  postSum_GenTime <- readRDS("PosteriorSummaries_GenTime_byArea.rds")
  postSum_VarDecomp <- readRDS("PosteriorSummaries_VarDecomp.rds")
  
  sheetNames <- c("HyperParameters", "Recruitment", "Survival", "RodentEffect", "PopDensity", "PopGrowthRate", "Detection")
  sheetMapping <- c("hParams.sum", "rRep.sum", "pSurv.sum", "betaR.sum", "popDens.sum", "lambda.sum", "detect.sum")
  
  
  sheetData <- list()
  
  for(i in 1:length(sheetNames)){
    
    idx <- which(names(postSum_main) == sheetMapping[i])
    
    sheetData[[i]] <-  postSum_main[[idx]]
    names(sheetData)[i] <- sheetNames[i]
  }
  
  openxlsx::write.xlsx(sheetData, file = "PosteriorSummaries_byArea.xlsx")
  
  
  readr::write_excel_csv(postSum_GenTime, "PosteriorSummaries_GenTime_byArea.csv")
  readr::write_excel_csv(postSum_VarDecomp, "PosteriorSummaries_VarDecomp.csv")
}
