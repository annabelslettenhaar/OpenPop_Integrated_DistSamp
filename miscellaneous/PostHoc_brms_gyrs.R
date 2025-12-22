
library(tidyverse)
library(sf)
library(terra)
library(parallel)
library(coda)

## Set localities/areas and time period of interest
areas <- c("Hardangervidda", 
           "Dovrefjell", 
           "Børgefjell")
minYear <- 1991
maxYear <- 2020

## Source all functions in "R_weather" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R_weather')

# Read in IDSM output
IDSM.out <- readRDS("/cloud/project/rypeIDSM_dHN_gyrData_23-10_mediumrun_fullloop_onlyOcc_RE.rds")


##########################################################################

## Read in and process weather data
d_weather <- wrangleData_GyrWeather(areas = areas,
                                    minYear = minYear, 
                                    maxYear = maxYear,
                                    byArea = FALSE)

# Define year range 
years <- seq_along(minYear:maxYear)

d_weather_long <- list()

# Loop through each element in the weather list
for (i in seq_along(d_weather)) {
  metric_name <- names(d_weather)[i]
  mat <- d_weather[[i]]$data
  
  # Convert to data frame
  df <- as.data.frame(mat)
  colnames(df) <- as.character(years)
  
  # Assign area index only if multiple rows (i.e., by_area = TRUE)
  if (nrow(df) > 1) {
    df$area <- seq_len(nrow(df))
  }
  
  # Reshape to long format
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = if ("area" %in% colnames(df)) colnames(df)[colnames(df) != "area"] else everything(),
      names_to = "year",
      values_to = "value"
    ) %>%
    mutate(year = as.integer(year)) %>%
    rename(!!metric_name := value)
  
  # Store in list
  d_weather_long[[metric_name]] <- df_long
}



## Extract posterior mean values and SD's for posthoc testing

epsOcc_summary <- extractParamSummary(modelOutput = IDSM.out,
                                      parameterName = "epsT.Occ",
                                      byArea = FALSE)

epsOcc_summary <- epsOcc_summary %>%
  mutate(meanOcc = mean,
         sdOcc = sd) %>%
  select(-c(mean, sd))

epsProd_summary <- extractParamSummary(modelOutput = IDSM.out,
                                       parameterName = "epsT.Prod",
                                       byArea = FALSE)

epsProd_summary <- epsProd_summary %>%
  mutate(meanProd = mean,
         sdProd = sd) %>%
  select(-c(mean, sd))

## Setup data for modelling

# Edit weather function to return means per year, not per year per area. 

residual_summary <- epsOcc_summary %>%
  left_join(epsProd_summary)

model_data <- residual_summary

for (df_name in names(d_weather_long)) {
  model_data <- model_data %>%
    left_join(d_weather_long[[df_name]], by = c("year"))
}


## Make correlation matrix

cor_matrix <- cor(model_data[, c("chicktemp", "apriltemp", "febtemp", "febsnow", "SD20", "precip", "longrain")], use = "complete.obs")

# Plot
library(ggcorrplot)
ggcorrplot(cor_matrix, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3, method = "square",
           colors = c("blue", "white", "red"))



## Modelling

library(brms)

# Occupancy
m1_occ <- brm(
  meanOcc | se(sdOcc) ~ febsnow + febtemp,
  data = model_data,
  family = gaussian()
)

plot(m1_occ)
summary(m3_temp)

pp_check(m1_occ)


# Productivity
m2_prod <- brm(
  meanProd | se(sdProd) ~ SD20 + chicktemp,
  data = model_data,
  family = gaussian()
)

P1 <- brm(
  meanProd | se(sdProd) ~ chicktemp + precip,
  data = model_data,
  family = gaussian()
)

P2 <- brm(
  meanProd | se(sdProd) ~ precip,
  data = model_data,
  family = gaussian()
)

P3 <- brm(
  meanProd | se(sdProd) ~ chicktemp + longrain,
  data = model_data,
  family = gaussian()
)

P4 <- brm(
  meanProd | se(sdProd) ~ chicktemp + precip,
  data = model_data,
  family = student()
)

plot(P4)
summary(m2_prod)

pp_check(P4, type = "dens_overlay")


ce <- conditional_effects(m2_temp)
plot(ce)

