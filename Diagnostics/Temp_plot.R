
## PLotting all the mean temperature data on a map

library(dplyr)
library(ggplot2)
library(sf)

plot_year <- 2015
startday <- 121
endday <- 153

temp_map <- d_temp_raw$df.temp_raw %>%
  mutate(
    Year = as.numeric(format(Date, "%Y")),
    JulianDay = as.numeric(format(Date, "%j"))
  ) %>%
  filter(
    Year == plot_year,
    JulianDay >= startday,
    JulianDay <= endday
  ) %>%
  group_by(grid_id) %>%
  summarise(
    temp_mean = mean(Value, na.rm = TRUE),
    .groups = "drop"
  )

grid_plot <- d_grid$grid_poly %>%
  left_join(temp_map, by = "grid_id")

ggplot(grid_plot) +
  geom_sf(aes(fill = temp_mean), colour = "grey70", linewidth = 0.1) +
  scale_fill_viridis_c(name = "Temp (°C)") +
  labs(
    title = paste0(
      "Mean temperature, JD ",
      startday, "-", endday,
      ", year ", plot_year
    )
  ) +
  theme_bw()