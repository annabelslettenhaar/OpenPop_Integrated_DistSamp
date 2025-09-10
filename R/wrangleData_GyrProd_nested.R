#' Prepare Gyrfalcon Productivity Data for Modeling
#'
#' Loads and reshapes chick productivity data from a local CSV file into
#' a 2 2D array with dimensions [area, year], one for the number of chicks produced
#' and one for the number of territories monitored.
#'
#' @param minYear Integer. The first year to include, used to index rows.
#' @param maxYear Integer. The last year to include, used to determine matrix size.
#'
#' @return A 3D numeric array [area, year, territory] of chick counts (or NA if missing).
#'         Dimension names are area names, years, and territory IDs.
#' @export


wrangleData_GyrProd_nested <- function(minYear, maxYear) {
  
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
    dplyr::select(Area, gyrArea, Year, TerritoryID, chicks, breeding_attempt) %>% 
    dplyr::filter(Year >= minYear, Year <= maxYear, breeding_attempt == 1) %>%
    dplyr::mutate(
      YearIdx   = Year - minYear + 1) %>%
    dplyr::filter(YearIdx >= 2)
  
  # Define separate arrays for chick counts, corresponding area and year
  chicksObs <- gyr_data$chicks
  chicksObs_area <- gyr_data$Area
  chicksObs_year <- gyr_data$YearIdx
  
  return(list(
    chicksObs   = chicksObs,
    chicksObs_area = chicksObs_area,
    chicksObs_year = chicksObs_year
  ))
}
