library(ggplot2)
library(cowplot)
library(dplyr)

# Get the data in the right format for plotting
samps <- as.matrix(onlyOcc)

# Convert to long format
long_df <- as.data.frame(samps) %>%
  tibble::rownames_to_column("iteration") %>%
  pivot_longer(-iteration, names_to = "param", values_to = "value") %>%
  filter(str_detect(param, "totDens_raw|GyrPressure_raw|S"))

# Extract variable type, area, and year
parsed <- str_match(long_df$param, "^([a-zA-Z_]+)\\[(\\d+),\\s*(\\d+)\\]$")
long_df$variable <- parsed[, 2]
long_df$area_index <- as.integer(parsed[, 3])
long_df$year_index <- as.integer(parsed[, 4])


# Map year/area index to actual values
long_df$year <- minYear + long_df$year_index
long_df <- long_df %>%
  mutate(year = minYear + year_index - 1,
         area = recode(area_index,
                       `1` = "Hardangervidda",
                       `2` = "Dovrefjell",
                       `3` = "Børgefjell"),
         value = if_else(variable == "totDens_raw", value * 1000000, value)) %>%
  filter(variable %in% c("totDens_raw","GyrPressure_raw","S"))

#Summarize posterior samples
summary_df <- long_df %>%
  group_by(variable, area, year) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(Variable = recode(variable,
                       "GyrPressure_raw" = "GyrPressure",
                       "S" = "Survival",
                       "totDens_raw" = "Ptarmigan density (birds/km^2)"))

# Plot with all three variables in one
p_prey <- ggplot(summary_df, aes(x = year, y = mean, color = Variable)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Variable),
              alpha = 0.2,
              color = NA) +  # removes the dark edge
  facet_grid(Variable ~ area, scales = "free_y") +
  scale_color_manual(values = c(
    "Ptarmigan density (birds/km^2)" = "#7EBB4C",     # green
    "GyrPressure" = "#E05D00", # orange
    "Survival" = "#CD348E"                # purple
  )) +
  scale_fill_manual(values = c(
    "Ptarmigan density (birds/km^2)" = "#7EBB4C",
    "GyrPressure" =  "#E05D00",
    "Survival" = "#CD348E"
  )) +
  theme_minimal(base_size = 14) + 
  theme(strip.text.y = element_blank(),
        legend.position = "bottom")
