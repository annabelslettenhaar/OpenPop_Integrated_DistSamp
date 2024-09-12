#' Summarises posteriors for average vital rates, population densities, and 
#' population growth rates per area
#'
#' @param mcmc.out an mcmc list containing posterior samples from a fitted model.
#' @param N_areas integer. Number of areas in the analysis. 
#' @param area_names character vector (length = N_areas) containing names of all
#' areas in the analysis. 
#' @param N_sites integer vector (length = N_areas) specifying the number of 
#' sites per area. 
#' @param min_years integer vector (length = N_areas) specifying the time index 
#' for the first year with data for each area. 
#' @param max_years integer vector (length = N_areas) specifying the time index 
#' for the first year with data for each area. 
#' @param minYear integer. First year included in analysis. 
#' @param maxYear integer. Last year included in analysis.
#' @param fitRodentCov logical. If TRUE, makes posterir summaries for rodent 
#' effect in addition to other parameters. 
#' @param save logical. If TRUE (default), saves the posterior summaries in an 
#' RDS file ()
#' 
#' @return a list containing 5 dataframes with summarised posteriors for several
#' parameters by area: average recruitment (rRep.sum), average survival 
#' (pSurv.sum), rodent effect slope (betaR.sum, = NA when fitRodentCov = FALSE), 
#' average population density (popDens.sum), and population growth rate 
#' (lambda.sum). 
#' @export
#'
#' @examples
#' 
summarisePost_areas <- function(mcmc.out, 
                                N_areas, area_names, N_sites, 
                                min_years, max_years, minYear, maxYear,
                                fitRodentCov, save = TRUE){
  
  
  ## Convert posterior samples to matrix format
  out.mat <- as.matrix(mcmc.out)
  
  ## Set year index ranges for plotting shared data collection period only
  minYearIdx_shared <- max(min_years)
  maxYearIdx_shared <- min(max_years)
  
  ## Determine area-specific year range
  area_yearIdxs <- (1:(maxYear-minYear+1))
  area_years <- area_yearIdxs + (minYear - 1)
  
  ## Calculate area- and year-specific population growth rates
  for(i in 1:N_areas){
    for(t in 1:length(area_years)){
      
      # Summarize annual average population densities (per area, averaged over lines)
      popDens_juv_t0 <- out.mat[, paste0("meanDens[",  i, ", 1, ", area_yearIdxs[t], "]")]
      popDens_ad_t0 <- out.mat[, paste0("meanDens[",  i, ", 2, ", area_yearIdxs[t], "]")]
      popDens_mean_t0 <- popDens_juv_t0 + popDens_ad_t0
      
      if(t < length(area_years)){
        popDens_juv_t1 <- out.mat[, paste0("meanDens[",  i, ", 1, ", area_yearIdxs[t+1], "]")]
        popDens_ad_t1 <- out.mat[, paste0("meanDens[",  i, ", 2, ", area_yearIdxs[t+1], "]")]
        popDens_mean_t1 <- popDens_juv_t1 + popDens_ad_t1
        
        # Calculate annual average population growth rate
        lambda_t0t1 <- popDens_mean_t1 / popDens_mean_t0
      }
      
      # Merge average density (summed age classes) back into posterior samples
      out.mat <- cbind(out.mat, popDens_mean_t0)
      colnames(out.mat)[dim(out.mat)[2]] <- paste0("meanDens[",i, ", ", area_yearIdxs[t], "]" )
      
      # Merge population growth rate back into posterior samples
      if(t < length(area_years)){
        out.mat <- cbind(out.mat, lambda_t0t1)
        colnames(out.mat)[dim(out.mat)[2]] <- paste0("lambda[",i, ", ", area_yearIdxs[t], "]" )
      }
    }
    
    ## Calculate average area-specific average densities and population growth rate over different time periods
    
    # Entire time period - average density
    lambda_mean_tot <- rowMeans(out.mat[, paste0("meanDens[", i, ", ", area_yearIdxs[1:(length(area_years)-1)], "]")])
    out.mat <- cbind(out.mat, lambda_mean_tot)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("densAvg_tot[", i, "]")
    
    # Overlapping time period (all areas) - average density
    lambda_mean_shared <- rowMeans(out.mat[, paste0("meanDens[", i, ", ", minYearIdx_shared:(maxYearIdx_shared-1), "]")])
    out.mat <- cbind(out.mat, lambda_mean_shared)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("densAvg_shared[", i, "]")
    
    # Entire time period - population growth rate
    lambda_mean_tot <- rowMeans(out.mat[, paste0("lambda[", i, ", ", area_yearIdxs[1:(length(area_years)-1)], "]")])
    out.mat <- cbind(out.mat, lambda_mean_tot)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("lambdaAvg_tot[", i, "]")
    
    # Overlapping time period (all areas) - population growth rate
    lambda_mean_shared <- rowMeans(out.mat[, paste0("lambda[", i, ", ", minYearIdx_shared:(maxYearIdx_shared-1), "]")])
    out.mat <- cbind(out.mat, lambda_mean_shared)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("lambdaAvg_shared[", i, "]")
  }
  
  ## Summarize posteriors for relevant parameters
  
  # Prepare matrices for storage of results
  rRep.sum <- pSurv.sum <- betaR.sum <- detect.sum <- data.frame()
  popDens.sum <- lambda.sum <- data.frame()
  
  for(i in 1:N_areas){
    
    # Summarize average reproductive rates
    rRep_name <- paste0("Mu.R[",  i, "]")
    rRep_add <- data.frame(Area = area_names[i],
                           Median = median(out.mat[, rRep_name]),
                           lCI = unname(quantile(out.mat[, rRep_name], probs = 0.025)),
                           uCI = unname(quantile(out.mat[, rRep_name], probs = 0.975)),
                           Mean = mean(out.mat[, rRep_name]),
                           SD = sd(out.mat[, rRep_name]),
                           CV = sd(out.mat[, rRep_name]) / mean(out.mat[, rRep_name])
    )
    rRep.sum <- rbind(rRep.sum, rRep_add)
    
    # Summarize average annual survival rates
    pSurv_name <- paste0("Mu.S[",  i, "]")
    pSurv_add <- data.frame(Area = area_names[i],
                            Median = median(out.mat[, pSurv_name]),
                            lCI = unname(quantile(out.mat[, pSurv_name], probs = 0.025)),
                            uCI = unname(quantile(out.mat[, pSurv_name], probs = 0.975)),
                            Mean = mean(out.mat[, pSurv_name]),
                            SD = sd(out.mat[, pSurv_name]),
                            CV = sd(out.mat[, pSurv_name])/mean(out.mat[, pSurv_name])
    )
    pSurv.sum <- rbind(pSurv.sum, pSurv_add)
    
    # Summarize rodent effects
    if(fitRodentCov){
      betaR_name <- paste0("betaR.R[",  i, "]")
      betaR_add <- data.frame(Area = area_names[i],
                              Median = median(out.mat[, betaR_name]),
                              lCI = unname(quantile(out.mat[, betaR_name], probs = 0.025)),
                              uCI = unname(quantile(out.mat[, betaR_name], probs = 0.975)),
                              Mean = mean(out.mat[, betaR_name]),
                              SD = sd(out.mat[, betaR_name]),
                              CV = sd(out.mat[, betaR_name]) / mean(out.mat[, betaR_name])
      )
      
      betaR.sum <- rbind(betaR.sum, betaR_add)
    } else{
      betaR.sum <- NA
    }
    
    # Summarize detection decay parameter
    detect_name <- paste0("mu.dd[",  i, "]")
    detect_add <- data.frame(Area = area_names[i],
                           Median = median(exp(out.mat[, detect_name])),
                           lCI = unname(quantile(exp(out.mat[, detect_name]), probs = 0.025)),
                           uCI = unname(quantile(exp(out.mat[, detect_name]), probs = 0.975)),
                           Mean = mean(exp(out.mat[, detect_name])),
                           SD = sd(exp(out.mat[, detect_name])),
                           CV = sd(exp(out.mat[, detect_name])) / mean(exp(out.mat[, detect_name])) 
    )
    detect.sum <- rbind(detect.sum, detect_add)
    
    # Summarize average population densities
    Dens_names <- c(paste0("densAvg_tot[",  i, "]"), paste0("densAvg_shared[",  i, "]"))
    Dens_tot_add <- data.frame(Area = area_names[i],
                               SummaryPeriod = paste0(minYear, "-", maxYear),
                               Median = median(out.mat[, Dens_names[1]]),
                               lCI = unname(quantile(out.mat[, Dens_names[1]], probs = 0.025)),
                               uCI = unname(quantile(out.mat[, Dens_names[1]], probs = 0.975)),
                               Mean = mean(out.mat[, Dens_names[1]]),
                               SD = sd(out.mat[, Dens_names[1]]),
                               CV = sd(out.mat[, Dens_names[1]]) / mean(out.mat[, Dens_names[1]])
    )
    Dens_shared_add <- data.frame(Area = area_names[i],
                                  SummaryPeriod = paste0(minYearIdx_shared + minYear - 1, "-", maxYearIdx_shared + minYear - 1),
                                  Median = median(out.mat[, Dens_names[2]]),
                                  lCI = unname(quantile(out.mat[, Dens_names[2]], probs = 0.025)),
                                  uCI = unname(quantile(out.mat[, Dens_names[2]], probs = 0.975)),
                                  Mean = mean(out.mat[, Dens_names[2]]),
                                  SD = sd(out.mat[, Dens_names[2]]),
                                  CV = sd(out.mat[, Dens_names[2]]) / mean(out.mat[, Dens_names[2]])
    )
    popDens.sum <- rbind(popDens.sum, Dens_tot_add, Dens_shared_add)
    
    # Summarize average population growth rates
    lambda_names <- c(paste0("lambdaAvg_tot[",  i, "]"), paste0("lambdaAvg_shared[",  i, "]"))
    lambda_tot_add <- data.frame(Area = area_names[i],
                                 SummaryPeriod = paste0(minYear, "-", maxYear),
                                 Median = median(out.mat[, lambda_names[1]]),
                                 lCI = unname(quantile(out.mat[, lambda_names[1]], probs = 0.025)),
                                 uCI = unname(quantile(out.mat[, lambda_names[1]], probs = 0.975)),
                                 Mean = mean(out.mat[, lambda_names[1]]),
                                 SD = sd(out.mat[, lambda_names[1]]),
                                 CV = sd(out.mat[, lambda_names[1]]) / mean(out.mat[, lambda_names[1]])
    )
    lambda_shared_add <- data.frame(Area = area_names[i],
                                    SummaryPeriod = paste0(minYearIdx_shared + minYear - 1, "-", maxYearIdx_shared + minYear - 1),
                                    Median = median(out.mat[, lambda_names[2]]),
                                    lCI = unname(quantile(out.mat[, lambda_names[2]], probs = 0.025)),
                                    uCI = unname(quantile(out.mat[, lambda_names[2]], probs = 0.975)),
                                    Mean = mean(out.mat[, lambda_names[2]]),
                                    SD = sd(out.mat[, lambda_names[2]]),
                                    CV = sd(out.mat[, lambda_names[2]]) / mean(out.mat[, lambda_names[2]])
    )
    lambda.sum <- rbind(lambda.sum, lambda_tot_add, lambda_shared_add)
  }
  
  ## Summarize hyper-parameters
  hParams.list <- c("h.Mu.R", "h.sigma.R", "sigmaT.R", "sigmaR.R",
                    "h.Mu.betaR.R", "h.sigma.betaR.R",
                    "h.Mu.S","h.sigma.S", "sigmaT.S", "sigmaR.S",
                    "h.mu.dd", "h.sigma.dd", "sigmaT.dd", "sigmaR.dd")
  hParams.sum <- data.frame()
  for(i in 1:length(hParams.list)){
    hParams_add <- data.frame(Parameter = hParams.list[i],
                              Median = median(out.mat[, hParams.list[i]]),
                              lCI = unname(quantile(out.mat[, hParams.list[i]], probs = 0.025)),
                              uCI = unname(quantile(out.mat[, hParams.list[i]], probs = 0.975)),
                              Mean = mean(out.mat[, hParams.list[i]]),
                              SD = sd(out.mat[, hParams.list[i]]),
                              CV = sd(out.mat[, hParams.list[i]]) / mean(out.mat[, hParams.list[i]])
    )
    hParams.sum <- rbind(hParams.sum, hParams_add)
  }
  
  ## Collect summary data into a list 
  PostSum.list <- list(rRep.sum = rRep.sum,
                       pSurv.sum = pSurv.sum,
                       betaR.sum = betaR.sum,
                       detect.sum = detect.sum,
                       popDens.sum = popDens.sum,
                       lambda.sum = lambda.sum,
                       hParams.sum = hParams.sum)
  
  ## Save (optional) and return
  if(save){
    saveRDS(PostSum.list, file = "PosteriorSummaries_byArea.rds")
  }

  return(PostSum.list)
}
  