library(tidyverse)
library(parallel)
library(coda)
library(ggcorrplot)

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


# Read in and process weather data

d_temp <- readRDS("data/weather/temp.rds") # Temperature during chick rearing

# d_temppre <- wrangleData_Temp(minYear = minYear,
#                               maxYear = maxYear,
#                               areas = areas,
#                               startday = 121,
#                               endday = 153)

d_temppre <- readRDS("data/weather/temppre.rds") # Temperature before and during incubation

# d_snow <- wrangleData_Snow(minYear = minYear,
#                            maxYear = maxYear,
#                            areas = areas)

d_snow <- readRDS("data/weather/snow.rds") # First 7-day snow free period

# d_snowdepth <- wrangleData_SnowDepth(minYear = minYear,
#                                      maxYear = maxYear,
#                                      areas = areas)

d_snowdepth <- readRDS("data/weather/snowdepth.rds") # Snow depth on the 20th of May

# d_tempJan <- wrangleData_Temp(minYear = minYear,
#                                  maxYear = maxYear,
#                                  areas = areas,
#                                  startday = 1,
#                                  endday = 31)

d_tempJan <- readRDS("data/weather/tempJan.rds") # Temperature in January

# d_snowperiod <- wrangleData_SnowPeriod(minYear = minYear,
#                                        maxYear = maxYear,
#                                        areas = areas)

d_snowperiod <- readRDS("data/weather/snowperiod.rds") # Total length of snow-free period

weather <- list(d_temp, d_temppre, d_tempJan, d_snow, d_snowdepth, d_snowperiod)

# Define year range 
years <- seq_along(minYear:maxYear)

# Initialize an empty list to store reshaped data frames
reshaped_list <- list()

# Define names for each metric (adjust to match the list structure)
metric_names <- c("TempPost", "TempPre", "TempJan", "SnowFree", "SD20", "SnowPeriod")

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



## Extract posterior mean values and SD's for posthoc testing

epsR_summary <- extractParamSummary(modelOutput = IDSM.out,
                                      parameterName = "epsR.R",
                                    byArea = TRUE)

epsR_summary <- epsR_summary %>%
  mutate(meanR = mean,
         sdR = sd) %>%
  select(-c(mean, sd))

epsS_summary <- extractParamSummary(modelOutput = IDSM.out,
                                    parameterName = "epsR.S",
                                    byArea = TRUE)

epsS_summary <- epsS_summary %>%
  mutate(meanS = mean,
         sdS = sd) %>%
  select(-c(mean, sd))

residual_summary <- epsR_summary %>%
  left_join(epsS_summary)

model_data <- residual_summary

for (df_name in names(reshaped_list)) {
  model_data <- model_data %>%
    left_join(reshaped_list[[df_name]], by = c("area", "year"))
}

model_data$area <- as.factor(model_data$area)


## Make correlation matrix

cor_matrix <- cor(model_data[, c("TempPost", "TempPre", "TempJan", "SnowFree", "SD20", "SnowPeriod")], use = "complete.obs")

# Plot
ggcorrplot(cor_matrix, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3, method = "square",
           colors = c("blue", "white", "red"))




## Modelling

library(brms)

# Recruitment

R1 <- brm(
  meanR | se(sdR) ~ TempPre + TempPost + SnowFree + (1 | area),
  data = model_data,
  family = gaussian()
)

R2 <- brm(
  meanR | se(sdR) ~ TempPre + TempPost + SD20 + (1 | area),
  data = model_data,
  family = gaussian()
)

R3 <- brm(
  meanR | se(sdR) ~ TempPre + TempPost + SnowFree + area,
  data = model_data,
  family = gaussian()
)

R4 <- brm(
  meanR | se(sdR) ~ TempPre + TempPost + area,
  data = model_data,
  family = gaussian()
)

R5 <- brm(
  meanR | se(sdR) ~ TempPre + SnowFree + area,
  data = model_data,
  family = gaussian()
)

R6 <- brm(
  meanR | se(sdR) ~ SnowFree + area,
  data = model_data,
  family = gaussian()
)

R7 <- brm(
  meanR | se(sdR) ~ TempPre + area,
  data = model_data,
  family = gaussian()
)

loo_R1 <- loo(R1)
loo_R2 <- loo(R2)
loo_R3 <- loo(R3)
loo_R4 <- loo(R4)
loo_R5 <- loo(R5)
loo_R6 <- loo(R6)
loo_R7 <- loo(R7)
loo_compare(loo_R1, loo_R2, loo_R3, loo_R4, loo_R6, loo_R5, loo_R7)

summary(R5)

library(emmeans)
pairs(emmeans(R3, ~ area))


plot(R7)
summary(R1)

pp_check(R7)

ce <- conditional_effects(R5)
plot(ce)

# Survival

S1 <- brm(
  meanS | se(sdS) ~ TempPre + TempJan + SnowPeriod + area,
  data = model_data,
  family = gaussian()
)

S2 <- brm(
  meanS | se(sdS) ~ TempPre + TempJan + SD20 + area,
  data = model_data,
  family = gaussian()
)

S3 <- brm(
  meanS | se(sdS) ~ TempJan + SnowPeriod + area,
  data = model_data,
  family = gaussian()
)

S4 <- brm(
  meanS | se(sdS) ~ TempPre + SnowPeriod + area,
  data = model_data,
  family = gaussian()
)

S5 <- brm(
  meanS | se(sdS) ~ TempPre + area,
  data = model_data,
  family = gaussian()
)

S6 <- brm(
  meanS | se(sdS) ~ TempJan + area,
  data = model_data,
  family = gaussian()
)

S7 <- brm(
  meanS | se(sdS) ~ SnowPeriod + area,
  data = model_data,
  family = gaussian()
)

loo_S1 <- loo(S1)
loo_S2 <- loo(S2)
loo_S3 <- loo(S3)
loo_S4 <- loo(S4)
loo_S5 <- loo(S5)
loo_S6 <- loo(S6)
loo_S7 <- loo(S7)
loo_compare(loo_S1, loo_S2, loo_S3, loo_S4, loo_S5, loo_S6, loo_S7)

plot(S4)
summary(S5)

pp_check(S7)


## Plot effects

library(dplyr)
library(tidyr)
library(stringr)

# Example for one model
samples_R5 <- posterior_samples(R5, pars = "^b_")
samples_R5 <- samples_R5 %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  mutate(model = "R5")

samples_S5 <- posterior_samples(S5, pars = "^b_")
samples_S5 <- samples_S5 %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  mutate(model = "S5")

samples_S6 <- posterior_samples(S6, pars = "^b_")
samples_S6 <- samples_S6 %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  mutate(model = "S6")

samples_S7 <- posterior_samples(S7, pars = "^b_")
samples_S7 <- samples_S7 %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  mutate(model = "S7")

effects_df <- bind_rows(samples_R5, samples_S5, samples_S6, samples_S7)  # Add all models

effects_df <- effects_df %>%
  mutate(
    type = str_replace(param, "b_", "")
  )

effects_df <- bind_rows(
  samples_R5 %>% mutate(model = "Recruitment 5"),
  samples_S5 %>% mutate(model = "Survival 5"),
  samples_S6 %>% mutate(model = "Survival 6"),
  samples_S7 %>% mutate(model = "Survival 7")
)

summary_df <- effects_df %>%
  group_by(model, type) %>%
  summarise(
    median = median(value),
    lci = quantile(value, 0.025),
    uci = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  filter(type %in% c("SnowFree", "SnowPeriod", "TempJan", "TempPre"))



library(ggplot2)

ggplot(summary_df, aes(x = median, y = type)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ model) +
  theme_minimal() +
  labs(x = "Posterior effect (median ± 95% CI)", y = "Covariate")
