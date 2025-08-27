# Extract total ptarmigan density per sqm and calculate means and sd per area
# For standardiziation in the model

extractTotalDensities <- function(mcmc.out, 
                                  N_areas, area_names, 
                                  minYear, maxYear) {
  
  out.mat <- as.matrix(mcmc.out)
  totDens_df <- data.frame()
  
  area_yearIdxs <- (1:(maxYear - minYear + 1))
  area_years <- area_yearIdxs + (minYear - 1)
  
  for(i in 1:N_areas){
    for(t in 1:length(area_years)){
      
      # Extract totDens directly from MCMC matrix
      tot_draws <- out.mat[, paste0("totDens_raw[", i, ", ", area_yearIdxs[t], "]")]
      
      # Summarize
      tot_add <- data.frame(
        Area   = i,
        Year   = area_years[t],
        Median = median(tot_draws),
        Mean   = mean(tot_draws),
        SD     = sd(tot_draws),
        lCI    = quantile(tot_draws, probs = 0.025),
        uCI    = quantile(tot_draws, probs = 0.975)
      )
      
      totDens_df <- rbind(totDens_df, tot_add)
    }
  }
  
  return(totDens_df)
}

# Example usage
totDens_df <- extractTotalDensities(
  mcmc.out   = IDSM.out, 
  N_areas    = input_data$nim.constant$N_areas, 
  area_names = input_data$nim.constant$area_names, 
  minYear    = minYear, 
  maxYear    = maxYear
)

head(totDens_df)
totDens_df$Area <- as.factor(totDens_df$Area)
ggplot(totDens_df, aes(x=Year, y=Median, colour=Area)) +
  geom_line()

totDens_mean <- totDens_df %>%
  group_by(Area) %>%
  summarise(areamean = mean(Mean, na.rm = TRUE),
            areasd = mean(SD, na.rm = TRUE))
  
totDens_mean
