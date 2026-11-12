## PLotting grid cell transect segment structure including individual observations

library(sf)
library(dplyr)
library(ggplot2)

focus_loc <- "EED285AC-0C02-466F-8E46-B3CE3D2F372C"

# Problematic observations
bad_obs <- d_grid$obs_segment_lookup %>%
  filter(
    locationID == focus_loc,
    is.na(segmentID)
  )

# All observations from this transect
obs_plot <- LT_data$d_obs_spatial %>%
  st_transform(25833) %>%
  filter(locationID == focus_loc) %>%
  left_join(
    bad_obs %>%
      select(eventID) %>%
      mutate(problem = TRUE),
    by = "eventID"
  ) %>%
  mutate(problem = ifelse(is.na(problem), FALSE, TRUE))

# Segments belonging to this transect
seg_plot <- d_grid$transect_segments %>%
  filter(locationID == focus_loc)

seg_plot %>%
  st_drop_geometry() %>%
  count(segmentID, grid_id)

# Relevant grid cells
grid_plot <- d_grid$grid_poly %>%
  filter(
    grid_id %in%
      unique(c(
        seg_plot$grid_id,
        bad_obs$grid_id
      ))
  )
ggplot() +
  
  geom_sf(
    data = grid_plot,
    fill = NA,
    colour = "grey50",
    linewidth = 0.7
  ) +
  
  # geom_sf(
  #   data = seg_plot,
  #   aes(colour = factor(segmentID)),
  #   linewidth = 1.5,
  #   show.legend = FALSE
  # ) +
  
  geom_sf(
    data = seg_plot,
    aes(colour = factor(grid_id)),
    linewidth = 2
  ) +
  
  geom_sf(
    data = obs_plot,
    aes(
      shape = problem,
      colour = problem
    ),
    size = 3
  ) +
  
  scale_colour_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "red"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "FALSE" = 16,
      "TRUE" = 17
    )
  ) +
  
  theme_bw() +
  
  ggtitle(
    paste(
      "Transect",
      focus_loc
    )
  )

st_geometry(d_grid$transect_segments)
names(d_grid$transect_segments)

d_grid$transect_segments %>%
  st_drop_geometry() %>%
  count(segmentID, grid_id) %>%
  count(segmentID) %>%
  filter(n > 1)

ggplot() +
  geom_sf(data = grid_plot,
          fill = NA,
          colour = "black") +
  geom_sf(
    data = seg_plot,
    aes(colour = factor(grid_id)),
    linewidth = 2
  ) +
  geom_sf_text(
    data = st_centroid(grid_plot),
    aes(label = grid_id)
  )

transect_plot <- LT_data$d_spatial %>%
  st_transform(25833) %>%
  filter(locationID == focus_loc)

ggplot() +
  geom_sf(data = grid_plot, fill = NA) +
  geom_sf(data = transect_plot,
          colour = "black",
          linewidth = 2) +
  geom_sf(data = seg_plot,
          aes(colour = factor(segmentID)),
          linewidth = 1)
