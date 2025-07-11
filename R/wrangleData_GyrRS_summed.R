#' Prepare Gyrfalcon Productivity Data for Modeling (Summed by Area-Year)
#'
#' Loads and reshapes chick productivity data from a local CSV file into
#' a 2D array with dimensions [area, year], where each value is the total number
#' of chicks produced across all territories in that area-year.
#'
#' @param minYear Integer. The first year to include, used to index rows.
#' @param maxYear Integer. The last year to include, used to determine matrix size.
#'
#' @return A 2D numeric array [area, year] of total chick counts (or NA if missing).
#'         Dimension names are area names and indexed years.
#' @export

wrangleData_GyrRS_summed <- function(minYear, maxYear) {
  
  if (!file.exists("data/Gyr_data.csv")) {
    stop("Data file (data/Gyr_data.csv) not found. This workflow requires this file.")
  }
  
  gyr_data_raw <- subset(read.csv("data/Gyr_data.csv"), select = -1)
  
  gyr_data <- gyr_data_raw %>%
    dplyr::select(Area, Year, TerritoryID, chicks) %>%
    dplyr::filter(Year >= minYear, Year <= maxYear) %>%
    dplyr::mutate(YearIdx = Year - minYear + 1)
  
  # Define dimensions
  areas <- sort(unique(gyr_data$Area))
  n_areas <- length(areas)
  n_years <- maxYear - minYear + 1
  
  # Initialize output array [area, year]
  array_out <- matrix(NA, nrow = n_areas, ncol = n_years,
                      dimnames = list(areas, paste0("Year", 1:n_years)))
  
  # Fill array with summed chick counts
  for (a in seq_along(areas)) {
    area_name <- areas[a]
    area_df <- gyr_data %>% dplyr::filter(Area == area_name)
    
    for (y in 1:n_years) {
      sum_val <- area_df %>%
        dplyr::filter(YearIdx == y) %>%
        dplyr::summarise(total_chicks = sum(chicks, na.rm = TRUE)) %>%
        dplyr::pull(total_chicks)
      
      # If all values were NA, keep it as NA
      if (all(is.na(area_df$chicks[area_df$YearIdx == y]))) {
        array_out[a, y] <- NA
      } else {
        array_out[a, y] <- sum_val
      }
    }
  }
  
  return(array_out)
}