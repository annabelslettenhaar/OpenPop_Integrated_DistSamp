#' Prepare Gyrfalcon Occupancy Data for Modeling
#'
#' Loads and reshapes breeding attempt (occupancy) data from a local CSV file 
#' into 2D arrays with dimensions [area, year]: 
#' one for the number of occupied territories, 
#' and one for the number of territories monitored.
#'
#' @param minYear Integer. The first year to include, used to index rows.
#' @param maxYear Integer. The last year to include, used to determine matrix size.
#'
#' @return A list with two 2D numeric arrays [area, year]:
#'   - Occ_count: number of occupied territories
#'   - n_monitored: number of monitored territories
#' @export


wrangleData_GyrOcc_agg <- function(minYear, maxYear) {
  
  if (!file.exists("data/Gyr_data.csv")) {
    stop("Data file (data/Gyr_data.csv) not found. This workflow requires this file.")
  }
  
  # Read and clean
  gyr_data_raw <- subset(read.csv("data/Gyr_data.csv"), select = -1)
  
  gyr_data_raw <- gyr_data_raw %>%
    mutate(gyrArea = recode(Area,
                            "1" = "Hardangervidda",
                            "2" = "Dovrefjell",
                            "3" = "Børgefjell"))
  
  gyr_data <- gyr_data_raw %>%
    dplyr::select(Area, gyrArea, Year, TerritoryID, breeding_attempt) %>% 
    dplyr::filter(Year >= minYear, Year <= maxYear) %>%
    dplyr::mutate(
      YearIdx   = Year - minYear + 1,
      monitored = ifelse(is.na(breeding_attempt), 0, 1)  # monitored if attempt was recorded
    )
  
  # Define dimensions
  areas   <- sort(unique(gyr_data$Area))
  n_areas <- length(areas)
  n_years <- maxYear - minYear + 1
  
  # Initialize arrays
  Occ_count    <- array(NA, dim = c(n_areas, n_years))
  n_monitored  <- array(NA, dim = c(n_areas, n_years))
  
  # Fill arrays by aggregating per area × year
  for (a in seq_along(areas)) {
    area_df <- gyr_data %>% filter(Area == areas[a])
    
    for (t in 1:n_years) {
      year_df <- area_df %>% filter(YearIdx == t)
      
      Occ_count[a, t]   <- sum(year_df$breeding_attempt, na.rm = TRUE)
      n_monitored[a, t] <- sum(year_df$monitored, na.rm = TRUE)
    }
  }
  
  # Drop dimension names (to match productivity function)
  dimnames(Occ_count)   <- NULL
  dimnames(n_monitored) <- NULL
  
  return(list(
    Occ_count   = Occ_count,
    n_monitored = n_monitored
  ))
}

