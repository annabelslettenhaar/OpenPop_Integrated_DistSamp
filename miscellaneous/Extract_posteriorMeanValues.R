
## Means per area, across all years

# Convert to matrix
samps <- as.matrix(IDSM.out)

# Grab only the totDens_raw variables
dens_cols <- grep("^epsR.R", colnames(samps), value = TRUE)

# Manually set the number of areas
N_areas <- 3   # <- change this to however many areas you have
areas <- 1:N_areas

# Container for results per area
results <- lapply(areas, function(a) {
  # subset columns for this area across all years
  cols_area <- grep(paste0("^epsR.R\\[", a, ","), dens_cols, value = TRUE)
  
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




## Means and sd per area, per year


# Convert to matrix and subset only epsR.R columns
samps <- as.matrix(IDSM.out)
epsR_cols <- grep("^epsR.R", colnames(samps), value = TRUE)

# Convert posterior samples to data frame
samps_df <- as.data.frame(samps[, epsR_cols])

# Add iteration index
samps_df$iteration <- seq_len(nrow(samps_df))

# Pivot longer
epsR_long <- samps_df %>%
  pivot_longer(
    cols = -iteration,
    names_to = "param",
    values_to = "value"
  )

# Extract area and year from column names
epsR_long <- epsR_long %>%
  mutate(
    param = as.character(param),
    area = as.integer(gsub("epsR.R\\[(\\d+),\\s*(\\d+)\\]", "\\1", param)),
    year = as.integer(gsub("epsR.R\\[(\\d+),\\s*(\\d+)\\]", "\\2", param))
  )

summary_epsR <- epsR_long %>%
  group_by(area, year) %>%
  summarise(
    epsR_mean = mean(value),
    epsR_sd = sd(value),
    .groups = "drop"
  )
