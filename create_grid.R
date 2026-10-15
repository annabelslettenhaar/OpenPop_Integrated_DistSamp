d_trans_grid <- LT_data$d_trans

if(areaAggregation){
  
  colnames(d_trans_grid)[
    which(colnames(d_trans_grid) == "verbatimLocality")
  ] <- "spatialUnit"
  
} else {
  
  colnames(d_trans_grid)[
    which(colnames(d_trans_grid) == "locality")
  ] <- "spatialUnit"
}


# Project transects
d_spatial_utm <- st_transform(LT_data$d_spatial, 25833)

# Create transect centroids
transect_centroids <- st_centroid(d_spatial_utm)

# Create grid
grid_all <- st_make_grid(
  d_spatial_utm,
  cellsize = 5000
) %>%
  st_sf()

# Keep only cells containing transects
cells_selected <- grid_all[
  lengths(st_intersects(grid_all, transect_centroids)) > 0,
]

# Assign unique grid IDs
cells_selected$grid_id <- 1:nrow(cells_selected)

# Grid centroids
grid_centroids <- st_centroid(cells_selected)

grid_xy <- st_coordinates(grid_centroids)

grid_df <- data.frame(
  grid_id = cells_selected$grid_id,
  x = grid_xy[,1],
  y = grid_xy[,2]
)

# Transect -> grid lookup
grid_area_lookup <- st_join(
  transect_centroids,
  cells_selected["grid_id"]
) %>%
  st_drop_geometry() %>%
  select(locationID, grid_id) %>%
  left_join(
    d_trans_grid %>%
      distinct(locationID, spatialUnit),
    by = "locationID"
  ) %>%
  distinct()

# Area assignment for grids
# grid_area <- grid_area_lookup %>%
#   group_by(grid_id) %>%
#   summarise(
#     spatialUnit = first(spatialUnit)
#   )
grid_area <- grid_area_lookup %>%
  distinct(grid_id, spatialUnit) # Like this each transect gets its own grid number

grid_df <- grid_df %>%
  left_join(grid_area, by = "grid_id")

