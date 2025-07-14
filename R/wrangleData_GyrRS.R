#' Prepare Gyrfalcon Productivity Data for Modeling
#'
#' Loads and reshapes chick productivity data from a local CSV file into
#' a 3D array with dimensions [area, year, territory].
#'
#' @param minYear Integer. The first year to include, used to index rows.
#' @param maxYear Integer. The last year to include, used to determine matrix size.
#'
#' @return A 3D numeric array [area, year, territory] of chick counts (or NA if missing).
#'         Dimension names are area names, years, and territory IDs.
#' @export

wrangleData_GyrRS <- function(minYear, maxYear) {
  
  if (!file.exists("data/Gyr_data.csv")) {
    stop("Data file (data/Gyr_data.csv) not found. This workflow requires this file.")
  }
  
  gyr_data_raw <- subset(read.csv("data/Gyr_data.csv"), select = -1)
  
  gyr_data_raw <- gyr_data_raw %>%
    mutate(gyrArea = recode(Area,
                            "1" = "Hardangervidda",
                            "2" = "Dovrefjell",
                            "3" = "Børgefjell"))
  
  gyr_data <- gyr_data_raw %>%
    dplyr::select(Area, gyrArea, Year, TerritoryID, chicks) %>%
    dplyr::filter(Year >= minYear, Year <= maxYear) %>%
    dplyr::filter(gyrArea %in% areas) %>%
    dplyr::mutate(YearIdx = Year - minYear + 1)
  
  # Define dimensions
  areas <- sort(unique(gyr_data$Area))
  territories <- sort(unique(gyr_data$TerritoryID))
  n_areas <- length(areas)
  n_territories <- length(territories)
  n_years <- maxYear - minYear + 1
  
  # Initialize array [area, year, territory]
  array_out <- array(NA, dim = c(n_areas, n_years, n_territories))
  
  # Fill array
  for (a in seq_along(areas)) {
    area_name <- areas[a]
    area_df <- gyr_data %>% dplyr::filter(Area == area_name)
    
    for (i in seq_len(nrow(area_df))) {
      row <- area_df[i, ]
      t_idx <- match(row$TerritoryID, territories)
      y_idx <- row$YearIdx
      array_out[a, y_idx, t_idx] <- row$chicks
    }
  }
  
  return(array_out)
}
