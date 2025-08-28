
wrangleData_GyrOcc_agg <- function(minYear, maxYear) {
  
  if (!file.exists("data/Gyr_data.csv")) {
    stop("Data file (data/Gyr_data.csv) not found. This workflow requires this file.")
  }
  
  gyr_data_raw <- subset(read.csv("data/Gyr_data.csv"), select = -1)
  
  gyr_data_raw <- gyr_data_raw %>%
    mutate(gyrArea = recode(Area,
                            "1" = "Hardangervidda",
                            "2" = "Dovrefjell",
                            "3" = "Børgefjell"))
  
  # Subset and clean data
  gyr_data <- gyr_data_raw %>%
    dplyr::select(Area, gyrArea, Year, TerritoryID, breeding_attempt) %>%
    dplyr::filter(Year >= minYear, Year <= maxYear)
  
  areas <- sort(unique(gyr_data$Area))
  n_areas <- length(areas)
  n_years <- maxYear - minYear + 1
  
  # Aggregate to area x year
  agg_df <- gyr_data %>%
    dplyr::group_by(Area, Year) %>%
    summarise(
      Occ_count = sum(breeding_attempt, na.rm = TRUE),         # number of territories occupied
      n_monitored = sum(!is.na(breeding_attempt)),             # number of territories monitored
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      YearIdx = Year - minYear + 1
    )
  
  # Optionally: reshape to wide array [area, year] if needed for nimble
  array_out <- array(NA, dim = c(n_areas, n_years, 2), 
                     dimnames = list(areas, 1:n_years, c("Occ_count", "n_monitored")))
  
  for (a in seq_along(areas)) {
    area_name <- areas[a]
    area_df <- agg_df %>% dplyr::filter(Area == area_name)
    
    for (i in seq_len(nrow(area_df))) {
      y_idx <- area_df$YearIdx[i]
      array_out[a, y_idx, "Occ_count"] <- area_df$Occ_count[i]
      array_out[a, y_idx, "n_monitored"] <- area_df$n_monitored[i]
    }
  }
  
  return(list(agg_df = agg_df, array_out = array_out))
}

