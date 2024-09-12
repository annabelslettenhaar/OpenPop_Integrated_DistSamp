#' Calculate posterior summaries of parameter correlations with density
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param N_areas integer. Number of areas included in analyses. 
#' @param area_names character vector containing area/location names.
#' @param N_sites matrix containing the number of sites per area/location.
#' @param min_years integer vector. Indices of first year with available data for each area/location. 
#' @param max_years integer vector. Indices of last year with available data for each area/location. 
#'
#' @return a vector of pdf plot names. The plots can be found in Plots/TimeSeries.
#' @export
#'
#' @examples

checkDD <- function(mcmc.out, 
                    N_areas, area_names, N_sites, 
                    min_years, max_years){
  
  ## Convert posterior samples to matrix format
  out.mat <- as.matrix(mcmc.out)

  ## Set up data frame for posterior summaries
  DD_corrs <- data.frame()
  
  ## Calculate posterior correlation coefficients 
  for(x in 1:N_areas){
    
    # Determine area-specific year range
    area_yearIdxs <- (min_years[x]:max_years[x])
    area_yearIdxs_surv <- (min_years[x]:(max_years[x]-1))
      
    # Set up vector for storing correlation coefficients
    corr_lambda <- corr_rep <- corr_surv <- rep(NA, nrow(out.mat))
    
    for(i in 1:nrow(out.mat)){
      
      # Extract parameter time series
      densT <- rep(NA, length(area_yearIdxs))

      for(t in 1:length(area_yearIdxs)){
        densT[t] <- out.mat[i, paste0("meanDens[",  x, ", 1, ", area_yearIdxs[t], "]")] + 
          out.mat[i, paste0("meanDens[",  x, ", 2, ", area_yearIdxs[t], "]")]
      }
      
      lambdaT <- c(densT[2:length(densT)]/densT[1:(length(densT)-1)], NA)
      repT <- unname(out.mat[i, paste0("R_year[",  x, ", ", area_yearIdxs, "]")])
      survT <- unname(out.mat[i, paste0("S[",  x, ", ", area_yearIdxs_surv, "]")])
      
      # Calculate correlation coefficients
      corr_lambda[i] <- cor.test(lambdaT, densT)$estimate
      corr_rep[i] <- cor.test(repT, densT)$estimate
      corr_surv[i] <- cor.test(survT, densT[1:(length(densT)-1)])$estimate
    }
    
    # Summarize posteriors in data frame
    corr_data <- rbind(quantile(corr_lambda, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                       quantile(corr_rep, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                       quantile(corr_surv, probs = c(0.5, 0.025, 0.25, 0.75, 0.975))) %>%
      as_tibble() %>%
      dplyr::mutate(Area = area_names[x],
                    Parameter = c("Lambda", "Recruitment", "Survival"),
                    N_years = length(area_yearIdxs),
                    .before = `50%`) %>%
      dplyr::mutate(Evidence = case_when(sign(`2.5%`) == sign(`97.5%`) ~ "**",
                                         sign(`25%`) == sign(`75%`) ~ "*",
                                         TRUE ~ "-"))
    DD_corrs <- rbind(DD_corrs, corr_data)
  }
  
  ## Save results as RDS
  saveRDS(DD_corrs, file = "DD_corrCoef.rds")
  write.csv(DD_corrs, "DD_corrCoef.csv", row.names = FALSE)
  
  ## Return filepaths
  filepaths <- c("DD_corrCoef.rds", "DD_corrCoef.csv")
}
