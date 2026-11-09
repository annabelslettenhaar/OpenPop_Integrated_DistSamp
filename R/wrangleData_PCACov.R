#' Extract PCA-based habitat covariates for model grid cells
#'
#' Aggregates PCA-derived habitat variables from a fine-resolution spatial layer
#' (e.g. 1 × 1 km habitat map) to the spatial grid used in the integrated
#' distance sampling model. The data originates from a PCA analysis by Bakkestuen
#' et al. 2008 (https://doi.org/10.1111/j.1365-2699.2008.01941.x). Maps used here 
#' are owned by NINA. PCA values are assigned to model grid cells based on
#' the location of each PCA cell centroid, and summary statistics are calculated
#' for every grid cell.
#'
#' The function returns mean and standard deviation of both PCA axes within each
#' model grid cell. Mean values are used directly as habitat covariates in
#' the model, whereas standard deviations provide a measure of habitat
#' heterogeneity within each grid cell.
#'
#' @param grid_poly An `sf` object containing the model grid polygons. Must
#'   include a column named `grid_id` with unique grid cell identifiers.
#' @param pca_map An `sf` object containing PCA-derived habitat data. Must
#'   contain columns named `PCA1` and `PCA2`.

wrangleData_PCACov <- function(grid_poly, pca_path = "data/PCA_maps") {
  
  library(sf)
  library(dplyr)
  
  ## Check whether PCA maps exists
  if (!dir.exists(pca_path)) {
    stop(
      paste0(
        "The PCA habitat map could not be found in '",
        pca_path,
        "'.\n\n",
        "The PCA habitat map is not distributed with this project ",
        "because it is based on data owned by The Norwegian Insitute ",
        "for Nature Research (NINA). ",
        "Please contact the data owners to obtain access and place ",
        "the shapefile components in the specified directory."
      )
    )
  }
  
  ## Locate shapefile
  shp_file <- list.files(
    pca_path,
    pattern = "\\.shp$",
    full.names = TRUE
  )
  
  if (length(shp_file) == 0) {
    stop(
      "No .shp file was found in the PCA map directory."
    )
  }
  
  ## Read map
  pca_map <- sf::st_read(shp_file, quiet = TRUE)
  

  
  ## Make sure grid dataframe and pca map have the same projection/CRS
  if (st_crs(grid_poly) != st_crs(pca_map)) {
    pca_map <- st_transform(pca_map, st_crs(grid_poly))
  }
  
  ## Get centroids of each grid cell in pca map
  pca_centroids <- st_centroid(pca_map)
  
  ## Overlay pca grid centroids with transect grid cells
  pca_grid <- st_join(pca_centroids, grid_poly["grid_id"], join = st_within, left = FALSE)
  
  ## Aggregate to grid cell level (means, and sd to get an idea of heterogeneity)
  grid_pca <- pca_grid %>%
    st_drop_geometry() %>%
    group_by(grid_id) %>%
    summarise(PCA1_mean = mean(PCA1, na.rm = TRUE),
              PCA2_mean = mean(PCA2, na.rm = TRUE),
              PCA1_sd   = sd(PCA1, na.rm = TRUE),
              PCA2_sd   = sd(PCA2, na.rm = TRUE),
              n_cells   = n(),
              .groups = "drop")%>%
    left_join(
      grid_poly |>
        st_drop_geometry() |>
        dplyr::select(grid_id, spatialUnit),
      by = "grid_id"
    ) %>%
    dplyr::relocate(spatialUnit, grid_id)
  
  ## Standardise means
  grid_pca <- grid_pca %>%
    mutate(PCA1 = as.numeric(scale(PCA1_mean)),
           PCA2 = as.numeric(scale(PCA2_mean)))
  
  ## Build x by g matrices for easy input in nimble data
  areas <- unique(grid_pca$spatialUnit)
  N_grid <- sapply(areas, function(a) sum(grid_pca$spatialUnit == a))
  max_grid <- max(N_grid)
  PCA1_mat <- matrix(NA, nrow = length(areas), ncol = max_grid, dimnames = list(areas, NULL))
  
  PCA2_mat <- matrix(NA, nrow = length(areas), ncol = max_grid, dimnames = list(areas, NULL))
  
  for (x in seq_along(areas)) { ## Filling the matrix
    tmp <- grid_pca %>%
      filter(spatialUnit == areas[x]) %>%
      arrange(grid_id)
    ng <- nrow(tmp)
    PCA1_mat[x, 1:ng] <- tmp$PCA1
    PCA2_mat[x, 1:ng] <- tmp$PCA2
  }
  
  return(list(d_pca = grid_pca,
              PCA1 = PCA1_mat,
              PCA2 = PCA2_mat
  ))
  
}


# grid_pca <- wrangleData_PCACov(
#   grid_poly = d_grid$grid_poly,
#   pca_map   = map_all
# )

# library(sf)
# 
# map_all <- st_read("data/PCA_maps/pca.shp")
# 
# cent <- st_centroid(map_all)
# 
# plot(
#   st_coordinates(cent),
#   col = hcl.colors(100)[cut(map_all$PCA1, 100)],
#   pch = 15,
#   cex = 0.5
# )
# 

