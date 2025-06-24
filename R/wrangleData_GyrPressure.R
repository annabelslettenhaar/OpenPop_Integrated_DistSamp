#' Prepare gyrfalcon occupancy data
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

wrangleData_GyrPressure <- function(localities = NULL, areas = NULL, areaAggregation, minYear, maxYear){
  
  ## Check if .csv file is available
  if(!file.exists("data/Gyr_data.csv")){
    stop("Data file (data/Gyr_data.csv) not found. The workflow currently requires this file as it does not yet support extraction of rodent observation data directly from GBIF/Living Norway.")
  }
  
  ## Load data from .csv file
  gyr_data_raw <- subset(read.csv("data/Gyr_data.csv"), select = -c(1))
  
  ## Rename areas to a more informative name & to match with ptarmigan datasets
  gyr_data_raw$Area <- as.factor(gyr_data_raw$Area)
  gyr_data_raw <- gyr_data_raw %>%
    mutate(gyrArea = recode(Area,
                            "1" = "Hardangervidda",
                            "2" = "Dovrefjell",
                            "3" = "Børgefjell"))
  
  ## Filter event data by either locality and year or area and year
  if(areaAggregation){
    gyr_data <- gyr_data_raw %>% 
      dplyr::filter(gyrArea %in% areas) %>%
      dplyr::filter(dplyr::between(Year, minYear, maxYear))
  }else{
    gyr_data <- gyr_data_raw %>% 
      dplyr::filter(TerritoryID %in% localities) %>%
      dplyr::filter(dplyr::between(Year, minYear, maxYear))
  }
  
  ## Double-check no duplicate transects remain
  duplicates <- gyr_data %>%
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
    colnames(gyr_data)[which(colnames(gyr_data) == "gyrArea")] <- "spatialUnit"
  }else{
    colnames(gyr_data)[which(colnames(gyr_data) == "TerritoryID")] <- "spatialUnit"
  }
  
  ## Calculate a measure of 'pressure' from the gyrfalcons that acts on ptarmigan survival
  gyr_data <- gyr_data %>%
    arrange(spatialUnit, TerritoryID, Year) %>%
    group_by(spatialUnit, TerritoryID) %>%
    mutate(breeding_attempt_minusone = lag(breeding_attempt)) # Occupancy previous year
  
  gyr_data <- gyr_data %>%
    mutate(pred_pressure = 0.5 *(2 * breeding_attempt_minusone) + 0.5 * (2 * breeding_attempt) + chicks)
  
  ## Summarise the pressure measure by spatial unit and year
  gyr_obs <- gyr_data %>% 
    dplyr::group_by(spatialUnit, Year) %>%
    dplyr::summarise(gyrPressure = sum(pred_pressure, na.rm = TRUE), .groups = "keep")
  
  ## Add year index
  gyr_obs$YearIdx <- gyr_obs$Year - minYear + 1
  
  ## Set up matrix for area-specific data
  gyrPressure <- matrix(NA, nrow = N_sUnits, ncol = length(minYear:maxYear))
  
  
  for(x in 1:N_sUnits){
    
    ## Subset data (specific area)
    if(!(sUnits[x] %in% gyr_obs$spatialUnit)){
      stop(paste0("Spatial unit ", sUnits[x], " (index ", x, ") is not in the data."))
    }
    
    gyr_obs_sub <- subset(gyr_obs, spatialUnit == sUnits[x])
    
    ## Extract year-specific data
    for(t in 1:length(minYear:maxYear)){
      
      if(t %in% gyr_obs_sub$YearIdx){
        gyrPressure[x, t] <- gyr_obs_sub$gyrPressure[which(gyr_obs_sub$YearIdx == t)]
      }
    }
  }
  
  ## Z-standardize covariate values
  meanCov <- mean(gyrPressure, na.rm = TRUE)
  sdCov <- sd(gyrPressure, na.rm = TRUE)
  #gyrOccAvg <- (gyrProdAvg - meanCov) / sdCov
  
  ## Return data
  return(list(gyrPressure = gyrPressure,
              meanCov = meanCov, 
              sdCov = sdCov))
  
}

# d_gyr_prod <- wrangleData_ProdGyr(#localities = localities,
#                             areas = areas,
#                             areaAggregation = areaAggregation,
#                             minYear = minYear, maxYear = maxYear)

