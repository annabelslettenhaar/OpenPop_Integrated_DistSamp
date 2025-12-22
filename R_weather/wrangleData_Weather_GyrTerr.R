#' Prepare weather data as covariates
#'
#' @param localities string or vector of strings. Names of localities to extract
#' data for. Either localities or areas must be provided. 
#' @param areas string or vector of strings. Names of areas to extract
#' data for. Either localities or areas must be provided.
#' @param areaAggregation logical. If TRUE, areas are used as smallest spatial unit. If FALSE, territories are used as smallest spatial unit.
#' @param minYear integer. Earliest year of data to extract.
#' @param maxYear integer. Latest year of data to extract.  
#'
#' @return a matrix containing the average number of transects with rodent observations per area and year.
#' @export
#'
#' @examples

# # For testing purposes
# localities <- listLocations()
# areas <- c("Hardangervidda", "Dovrefjell", "Børgefjell")
# minYear <- 1991
# maxYear <- 2020

wrangleData_Weather_GT <- function(localities = NULL, areas = NULL, areaAggregation, minYear, maxYear){
  
  ## Check if .csv file is available
  if(!file.exists("data/Gyr_data.csv")){
    stop("Data file (data/Gyr_data.csv) not found. The workflow currently requires this file as it does not yet support extraction of rodent observation data directly from GBIF/Living Norway.")
  }
  
  ## Load data from .csv file
  weather_data_raw <- subset(read.csv("data/gyrnestID_weather_ptar_cov.csv"), select = -c(1))
  
  ## Rename areas to a more informative name & to match with ptarmigan datasets
  weather_data_raw$Area <- as.factor(weather_data_raw$Area)
  weather_data_raw <- weather_data_raw %>%
    mutate(gyrArea = recode(Area,
                         "1" = "Hardangervidda",
                         "2" = "Dovrefjell",
                         "3" = "Børgefjell"))
  
  ## Filter event data by either locality and year or area and year
  if(areaAggregation){
    weather_data <- weather_data_raw %>% 
      dplyr::filter(gyrArea %in% areas) %>%
      dplyr::filter(dplyr::between(Year, minYear, maxYear))
  }else{
    weather_data <- weather_data_raw %>% 
      dplyr::filter(TerritoryID %in% localities) %>%
      dplyr::filter(dplyr::between(Year, minYear, maxYear))
  }
  
  ## Double-check no duplicate transects remain
  duplicates <- weather_data %>%
    dplyr::group_by(TerritoryID, gyrArea, Year) %>%
    dplyr::summarise(observations = dplyr::n(), .groups = 'keep') %>%
    dplyr::filter(observations > 1)
  
  if(nrow(duplicates) > 0){
    stop("There are duplicate transects (> 1 transect in same location per year).")
  }
  
  ## Assignment of spatial units
  if(areaAggregation){
    sUnits <- areas
  }else{
    sUnits <- localities
  }
  N_sUnits <- length(sUnits)
  
  ## Rename appropriate column in line transect data to reflect level of spatial aggregation
  if(areaAggregation){
    colnames(weather_data)[which(colnames(weather_data) == "gyrArea")] <- "spatialUnit"
  }else{
    colnames(weather_data)[which(colnames(weather_data) == "TerritoryID")] <- "spatialUnit"
  }
  
  ## Function to create standardized lists of the raw weather values, mean and sd for each variable
  setupWeatherData <- function(df, var, sUnits, N_sUnits, minYear, maxYear) {
    # Summarize mean value per spatialUnit and year
    df_sum <- df %>%
      group_by(spatialUnit, Year) %>%
      summarise(val = mean(.data[[var]], na.rm = TRUE), .groups = "keep") %>%
      mutate(YearIdx = Year - minYear + 1)
    
    # Create data matrix
    mat <- matrix(NA, nrow = N_sUnits, ncol = length(minYear:maxYear))
    for (x in seq_len(N_sUnits)) {
      unit <- sUnits[x]
      df_sub <- df_sum[df_sum$spatialUnit == unit, ]
      for (t in seq_len(ncol(mat))) {
        if (t %in% df_sub$YearIdx) {
          mat[x, t] <- df_sub$val[df_sub$YearIdx == t]
        }
      }
    }
    
    # Compute mean and sd for standardization
    cov_mean <- mean(mat, na.rm = TRUE)
    cov_sd   <- sd(mat, na.rm = TRUE)
    
    # Standardize
    mat_std <- (mat - cov_mean) / cov_sd
    
    list(
      data = mat_std,
      mean = cov_mean,
      sd = cov_sd
    )
  }
  
  # Build each covariate list
  d_SD <- setupWeatherData(weather_data, "snow_prebrood", sUnits, N_sUnits, minYear, maxYear)
  d_spring <- setupWeatherData(weather_data, "snowfree_jd", sUnits, N_sUnits, minYear, maxYear)
  d_winter <- setupWeatherData(weather_data, "snowstart_jd", sUnits, N_sUnits, minYear, maxYear)
  d_temp <- setupWeatherData(weather_data, "temp_brood", sUnits, N_sUnits, minYear, maxYear)
  
  return(list(
    d_SD     = d_SD,
    d_spring = d_spring,
    d_winter = d_winter,
    d_temp = d_temp
  ))
}

# weather_data <- wrangleData_Weather(areas = areas,
#                                     areaAggregation = areaAggregation,
#                                     minYear = minYear, maxYear = maxYear)
