# After running a testround / full round of the model, use the code below to extract population density data

extractPopulationDensities <- function(mcmc.out, 
                                       N_areas, area_names, 
                                       min_years, max_years, 
                                       minYear, maxYear) {
  
  out.mat <- as.matrix(mcmc.out)
  popDens <- data.frame()
  
  for(i in 1:N_areas){
    
    area_yearIdxs <- (1:(maxYear - minYear + 1))
    area_years <- area_yearIdxs + (minYear - 1)
    
    for(t in 1:length(area_years)){
      
      # Extract juvenile and adult densities
      popDens_juv <- out.mat[, paste0("meanDens[", i, ", 1, ", area_yearIdxs[t], "]")]
      popDens_ad  <- out.mat[, paste0("meanDens[", i, ", 2, ", area_yearIdxs[t], "]")]
      popDens_mean <- popDens_juv + popDens_ad
      
      # Summarize
      popDens_add <- data.frame(
        Area = area_names[i],
        Year = area_years[t],
        Median = median(popDens_mean),
        lCI = quantile(popDens_mean, probs = 0.025),
        uCI = quantile(popDens_mean, probs = 0.975)
      )
      
      popDens <- rbind(popDens, popDens_add)
    }
  }
  
  return(popDens)
}


popDens_df <- extractPopulationDensities(
  mcmc.out = density.out.tidy, 
  N_areas = input_data$nim.constant$N_areas, 
  area_names = input_data$nim.constant$area_names, 
  min_years = input_data$nim.constant$min_years, 
  max_years = input_data$nim.constant$max_years, 
  minYear = minYear, maxYear = maxYear)

head(popDens_df)

write.csv(popDens_df, "data/ptar_density_placeholder.csv")
