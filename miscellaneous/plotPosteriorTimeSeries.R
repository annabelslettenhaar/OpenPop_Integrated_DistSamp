
## Plot and summarise posterior samples

library(coda)
library(tidyverse)
library(cowplot)

posteriorsamples <- as.matrix(IDSM.out)

# Convert to long format
long_df <- as.data.frame(posteriorsamples) %>%
  tibble::rownames_to_column("iteration") %>%
  pivot_longer(-iteration, names_to = "param", values_to = "value") %>%
  filter(str_detect(param, "totDens_raw|GyrPressure_raw"))

long_df <- as.data.frame(posteriorsamples) %>%
  tibble::rownames_to_column("iteration") %>%
  pivot_longer(-iteration, names_to = "param", values_to = "value") %>%
  filter(str_detect(param, "probOcc|terrProd"))

long_df <- as.data.frame(posteriorsamples) %>%
  tibble::rownames_to_column("iteration") %>%
  pivot_longer(-iteration, names_to = "param", values_to = "value") %>%
  filter(str_detect(param, "^S\\[|^R_year\\["))



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
         value = if_else(variable == "totDens_raw", value * 1000000, value)
  ) 

#Summarize posterior samples
summary_df <- long_df %>%
  group_by(variable, area, year) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  )

# Plot

# Define a scaling factor to align the second variable
hv <- summary_df %>% filter(area == "Hardangervidda")

scale_factor_hv <- max(hv$mean[hv$variable == "terrProd"]) /
  max(hv$mean[hv$variable == "probOcc"])

# Define consistent labels
hv <- hv %>%
  mutate(variable_label = case_when(
    variable == "terrProd" ~ "Productivity per occupied territory",
    variable == "probOcc" ~ "Probability of occupancy"
  ))

hv_plot <- ggplot() +
  # Prey density
  geom_line(data = hv %>% filter(variable == "terrProd"),
            aes(x = year, y = mean, color = variable_label)) +
  geom_ribbon(data = hv %>% filter(variable == "terrProd"),
              aes(x = year, ymin = lower, ymax = upper, fill = variable_label), alpha = 0.2) +
  # Predator occupancy
  geom_line(data = hv %>% filter(variable == "probOcc"),
            aes(x = year, y = mean * scale_factor_hv, color = variable_label)) +
  geom_ribbon(data = hv %>% filter(variable == "probOcc"),
              aes(x = year, ymin = lower * scale_factor_hv, ymax = upper * scale_factor_hv, fill = variable_label),
              alpha = 0.2) +
  
  scale_y_continuous(
    name = "Productivity per occupied territory",
    sec.axis = sec_axis(~ . / scale_factor_hv, name = "Probability of occupancy")
  ) +
  scale_color_manual(name = "Variable",
                     values = c("Productivity per occupied territory" = "darkgreen", "Probability of occupancy" = "#d95f02")) +
  scale_fill_manual(name = "Variable",
                    values = c("Productivity per occupied territory" = "darkgreen", "Probability of occupancy" = "#d95f02")) +
  labs(title = "Posterior estimates for Hardangervidda", 
       x = "Year") +
  theme_minimal()



hv_plot <- ggplot() +
  # Prey density
  geom_line(data = hv %>% filter(variable == "terrProd"),
            aes(x = year, y = mean, color = "Productivity")) +
  geom_ribbon(data = hv %>% filter(variable == "terrProd"),
              aes(x = year, ymin = lower, ymax = upper, fill = "Productivity"), alpha = 0.2) +
  # Predator occupancy
  geom_line(data = hv %>% filter(variable == "probOcc"),
            aes(x = year, y = mean * scale_factor_hv, color = "Probability of occupancy")) +
  geom_ribbon(data = hv %>% filter(variable == "probOcc"),
              aes(x = year, ymin = lower * scale_factor_hv, ymax = upper * scale_factor_hv, fill = "Probability of occupancy"),
              alpha = 0.2) +
  
  scale_y_continuous(
    name = "Prey Density",
    sec.axis = sec_axis(~ . / scale_factor_hv, name = "Predator Occupancy")
  ) +
  scale_color_manual(values = c("Prey Density" = "darkgreen", "Predator Occupancy" = "#d95f02")) +
  scale_fill_manual(values = c("Prey Density" = "darkgreen", "Predator Occupancy" = "#d95f02")) +
  labs(title = "Posterior estimates for Hardangervidda", 
       x = "Year", color = "Variable", fill = "Variable") +
  theme_minimal()


df <- summary_df %>% filter(area == "Dovrefjell")
scale_factor_df <- max(df$mean[df$variable == "totDens_raw"]) /
  max(df$mean[df$variable == "GyrPressure_raw"])

df_plot <- ggplot() +
  # Prey density
  geom_line(data = df %>% filter(variable == "totDens_raw"),
            aes(x = year, y = mean, color = "Prey Density")) +
  geom_ribbon(data = df %>% filter(variable == "totDens_raw"),
              aes(x = year, ymin = lower, ymax = upper, fill = "Prey Density"), alpha = 0.2) +
  # Predator occupancy
  geom_line(data = df %>% filter(variable == "GyrPressure_raw"),
            aes(x = year, y = mean * scale_factor_df, color = "Predator Occupancy")) +
  geom_ribbon(data = df %>% filter(variable == "GyrPressure_raw"),
              aes(x = year, ymin = lower * scale_factor_df, ymax = upper * scale_factor_df, fill = "Predator Occupancy"),
              alpha = 0.2) +
  
  scale_y_continuous(
    name = "Prey Density",
    sec.axis = sec_axis(~ . / scale_factor_df, name = "Predator Occupancy")
  ) +
  scale_color_manual(values = c("Prey Density" = "darkgreen", "Predator Occupancy" = "#d95f02")) +
  scale_fill_manual(values = c("Prey Density" = "darkgreen", "Predator Occupancy" = "#d95f02")) +
  labs(title = "Posterior estimates for Dovrefjell", 
       x = "Year", color = "Variable", fill = "Variable") +
  theme_minimal()


bf <- summary_df %>% filter(area == "Børgefjell")
scale_factor_bf <- max(hv$mean[bf$variable == "totDens_raw"]) /
  max(hv$mean[bf$variable == "GyrPressure_raw"])

bf_plot <- ggplot() +
  # Prey density
  geom_line(data = bf %>% filter(variable == "totDens_raw"),
            aes(x = year, y = mean, color = "Prey Density")) +
  geom_ribbon(data = bf %>% filter(variable == "totDens_raw"),
              aes(x = year, ymin = lower, ymax = upper, fill = "Prey Density"), alpha = 0.2) +
  # Predator occupancy
  geom_line(data = bf %>% filter(variable == "GyrPressure_raw"),
            aes(x = year, y = mean * scale_factor_bf, color = "Predator Occupancy")) +
  geom_ribbon(data = bf %>% filter(variable == "GyrPressure_raw"),
              aes(x = year, ymin = lower * scale_factor_bf, ymax = upper * scale_factor_bf, fill = "Predator Occupancy"),
              alpha = 0.2) +
  
  scale_y_continuous(
    name = "Prey Density",
    sec.axis = sec_axis(~ . / scale_factor_bf, name = "Predator Occupancy")
  ) +
  scale_color_manual(values = c("Prey Density" = "darkgreen", "Predator Occupancy" = "#d95f02")) +
  scale_fill_manual(values = c("Prey Density" = "darkgreen", "Predator Occupancy" = "#d95f02")) +
  labs(title = "Posterior estimates for Børgefjell", 
       x = "Year", color = "Variable", fill = "Variable") +
  theme_minimal()


timeseries_allareas <- plot_grid(
  hv_plot, df_plot, bf_plot,
  ncol = 1,         # Stack vertically
  align = "v",      # Align vertically
  axis = "lr"       # Align left and right axes
)
