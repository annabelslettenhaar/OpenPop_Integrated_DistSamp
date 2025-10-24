
# Step 1: Convert matrix to data frame
temp_df <- as.data.frame(d_temp$data)
snow_df <- as.data.frame(d_snow$data)

# Step 2: Assign column names as actual years (e.g., 1996 to 2025)
# Adjust this range to match your actual years
years <- 1:ncol(temp_df)
years <- 1:ncol(snow_df) # or use actual years like 1996:2025
colnames(temp_df) <- as.character(years)
colnames(snow_df) <- as.character(years)

# Step 3: Add area identifiers
temp_df$area <- seq_len(nrow(temp_df))
snow_df$area <- seq_len(nrow(snow_df))

# Step 4: Reshape to long format
temp_long <- temp_df %>%
  pivot_longer(
    cols = -area,
    names_to = "year",
    values_to = "TempMean"
  ) %>%
  mutate(year = as.integer(year))

snow_long <- snow_df %>%
  pivot_longer(
    cols = -area,
    names_to = "year",
    values_to = "SnowFree"
  ) %>%
  mutate(year = as.integer(year))


data_for_model <- summary_epsR %>%
  left_join(temp_long, by = c("area", "year")) %>%
  left_join(snow_long, by = c("area", "year"))

plot(data_for_model$TempMean, data_for_model$SnowFree)

library(brms)

m3_snow <- brm(
  epsR_mean | se(epsR_sd) ~ SnowFree + (1 | area),
  data = data_for_model,
  family = gaussian()
)

plot(m3_snow)
summary(m3_tempsnow)

pp_check(m2_temp)


ce <- conditional_effects(m2_temp)
plot(ce)

