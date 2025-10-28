
# Read in IDSM output
IDSM.out <- readRDS("/cloud/project/rypeIDSM_dHN_gyrData_23-10_mediumrun_fullloop_onlyOcc_RE.rds")


##########################################################################

## RECRUITMENT


# Read in weather data
d_temp <- readRDS("data/weather/temp.rds")

# d_temppre <- wrangleData_Temp(minYear = minYear,
#                               maxYear = maxYear,
#                               areas = areas,
#                               startday = 121,
#                               endday = 153)
# saveRDS(d_temppre, "data/weather/temppre.rds")
d_temppre <- readRDS("data/weather/temppre.rds")

# d_snow <- wrangleData_Snow(minYear = minYear,
#                            maxYear = maxYear,
#                            areas = areas)
# saveRDS(d_snow, "data/weather/snow.rds")
d_snow <- readRDS("data/weather/snow.rds")

# d_snowdepth <- wrangleData_SnowDepth(minYear = minYear,
#                                      maxYear = maxYear,
#                                      areas = areas)
# saveRDS(d_snowdepth, "data/weather/snowdepth.rds")
d_snowdepth <- readRDS("data/weather/snowdepth.rds")

weather <- list(d_snow, d_snowdepth, d_temp, d_temppre)


# Define year range 
years <- 1:30

# Initialize an empty list to store reshaped data frames
reshaped_list <- list()

# Define names for each metric (adjust to match your list structure)
metric_names <- c("SnowFree", "SD20", "TempPost", "TempPre")

# Loop through each element in the weather list
for (i in seq_along(weather)) {
  # Choose which matrix to use: 'data' or 'standardized'
  if (!is.null(weather[[i]]$standardized)) {
    mat <- weather[[i]]$standardized
  } else {
    mat <- weather[[i]]$data
  }
  
  # Convert to data frame
  df <- as.data.frame(mat)
  colnames(df) <- as.character(years)
  df$area <- rownames(mat) %||% seq_len(nrow(df))  # fallback if no rownames
  
  # Reshape to long format
  df_long <- df %>%
    pivot_longer(
      cols = -area,
      names_to = "year",
      values_to = metric_names[i]
    ) %>%
    mutate(year = as.integer(year))
  
  # Store in list
  reshaped_list[[metric_names[i]]] <- df_long
}




# Extract posterior mean values and SD's for posthoc testing
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


data_for_model <- summary_epsR

for (df_name in names(reshaped_list)) {
  data_for_model <- data_for_model %>%
    left_join(reshaped_list[[df_name]], by = c("area", "year"))
}


data_for_model <- summary_epsR %>%
  left_join(temp_long, by = c("area", "year")) %>%
  left_join(snow_long, by = c("area", "year"))

plot(data_for_model$SnowFree, data_for_model$SD20)

## Make correlation matrix

cor_matrix <- cor(data_for_model[, c("SnowFree", "SD20", "TempPost", "TempPre")], use = "complete.obs")

# Plot
library(ggcorrplot)
ggcorrplot(cor_matrix, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3, method = "square",
           colors = c("blue", "white", "red"))



library(brms)

m3_temp <- brm(
  epsR_mean | se(epsR_sd) ~ TempPost + TempPre + (1 | area),
  data = data_for_model,
  family = gaussian()
)

plot(m3_temp)
summary(m3_temp)

pp_check(m2_temp)


ce <- conditional_effects(m2_temp)
plot(ce)


##########################################################################

## SURVIVAL




