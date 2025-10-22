
library(coda)

# Convert to matrix
samps <- as.matrix(IDSM.out)

# Grab only the totDens_raw variables
dens_cols <- grep("^GyrPressure", colnames(samps), value = TRUE)

# Manually set the number of areas
N_areas <- 3   # <- change this to however many areas you have
areas <- 1:N_areas

# Container for results per area
results <- lapply(areas, function(a) {
  # subset columns for this area across all years
  cols_area <- grep(paste0("^GyrPressure\\[", a, ","), dens_cols, value = TRUE)
  
  # compute mean and sd across years per posterior draw
  mean_across <- rowMeans(samps[, cols_area, drop = FALSE])
  sd_across   <- apply(samps[, cols_area, drop = FALSE], 1, sd)
  
  data.frame(
    area = a,
    mean_across_years = mean_across,
    sd_across_years   = sd_across
  )
})

# Combine all areas into one data frame
posterior_summary <- do.call(rbind, results)

par(mfrow = c(2, N_areas))  # 2 rows: means & SDs, N_areas columns

for(a in 1:N_areas){
  subset_a <- posterior_summary[posterior_summary$area == a, ]
  
  # histogram of mean across years
  hist(subset_a$mean_across_years, breaks = 40,
       main = paste("Mean (Area", a, ")"),
       xlab = "Mean totDens_raw")
  
  # histogram of SD across years
  hist(subset_a$sd_across_years, breaks = 40,
       main = paste("SD (Area", a, ")"),
       xlab = "SD totDens_raw")
}

median_means_per_area <- aggregate(mean_across_years ~ area, data = posterior_summary, median)
median_sd_per_area <- aggregate(sd_across_years ~ area, data = posterior_summary, median)

