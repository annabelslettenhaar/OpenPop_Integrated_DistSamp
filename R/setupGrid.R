#' Set up a spatial grid for transect-based analyses
#'
#' @param d_trans tibble containing information on transects (events) and their location names
#' @param d_spatial tibble containing transect geometries
#' @param areaAggregation logical. If TRUE, areas are used as smallest spatial unit. 
#' If FALSE, locations (within areas) are used as smallest spatial unit.
#' @param gridOnlyTransects Logical. If TRUE, only grid cells containing
#' at least one transect centroid are retained. If FALSE, all grid cells
#' covering the study area are retained.
#' @param grid_cell_size Numeric. Grid cell size in metres.
#'
#'
#' @return A list of dataframes: one containing centroid coordinates and their corresponding
#' spatialUnit name (grid_df), and one containing information that links transects to grid
#' cells and their spatial Units (grid_area_lookup)

setupGrid <- function(d_trans, d_spatial, 
                      areaAggregation, gridOnlyTransects,
                      grid_cell_size = 5000
                      ) {
  
  # Transect data
  d_trans_grid <- d_trans
  
  if (areaAggregation) {
    colnames(d_trans_grid)[colnames(d_trans_grid) == "verbatimLocality"] <- "spatialUnit"
  } else {
    colnames(d_trans_grid)[colnames(d_trans_grid) == "locality"] <- "spatialUnit"
  }
  
  # Project transects
  d_spatial_utm <- sf::st_transform(d_spatial, 25833)
  
  # Create transect centroids
  transect_centroids <- sf::st_centroid(d_spatial_utm)
  
  # Create grid
  grid_all <- sf::st_make_grid(d_spatial_utm, cellsize = grid_cell_size) %>%
    sf::st_sf()
  
  # Select grid cells
  if (gridOnlyTransects) {
    # Keep only cells containing transect centroids
    cells_selected <- grid_all[lengths(st_intersects(grid_all,transect_centroids)) > 0,
    ]
  } else {
    # Keep entire grid
    cells_selected <- grid_all
  }
  
  # Assign unique grid IDs
  cells_selected$grid_id <- seq_len(nrow(cells_selected))
  
  # Grid centroids
  grid_centroids <- sf::st_centroid(cells_selected)
  
  grid_xy <- sf::st_coordinates(grid_centroids)
  
  grid_df <- data.frame(
    grid_id = cells_selected$grid_id,
    x = grid_xy[, 1],
    y = grid_xy[, 2]
  )
  
  # Transect -> grid lookup
  grid_area_lookup <- sf::st_join(transect_centroids, cells_selected["grid_id"]) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(locationID, grid_id) %>%
    dplyr::left_join(d_trans_grid %>%
        dplyr::distinct(locationID, spatialUnit), by = "locationID") %>%
    dplyr::distinct()
  
  # Area assignment for grids
  grid_area <- grid_area_lookup %>%
    dplyr::distinct(grid_id, spatialUnit)
  
  grid_df <- grid_df %>%
    dplyr::left_join(grid_area, by = "grid_id")
  
  cells_selected <- cells_selected %>%
    dplyr::left_join(grid_area, by = "grid_id")
  
  return(list(
      grid_df = grid_df,
      grid_area_lookup = grid_area_lookup,
      grid_poly = cells_selected
    )
  )
}


# d_grid <- setupGrid(d_trans = LT_data$d_trans, 
#                     d_spatial = LT_data$d_spatial,
#                     areaAggregation = areaAggregation,
#                     gridOnlyTransects = gridOnlyTransects,
#                     grid_cell_size = 5000)
