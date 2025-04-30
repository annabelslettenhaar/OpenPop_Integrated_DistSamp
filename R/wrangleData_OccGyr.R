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

## For testing purposes
#areaAggregation <- TRUE
#areas <- c(1, 2, 3)
#localities <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#minYear <- 1991
#maxYear <- 2020

wrangleData_OccGyr <- function(localities = NULL, areas = NULL, areaAggregation, minYear, maxYear){

  ## Check if .csv file is available
  if(!file.exists("data/Gyr_data.csv")){
    stop("Data file (data/Gyr_data.csv) not found. The workflow currently requires this file as it does not yet support extraction of rodent observation data directly from GBIF/Living Norway.")
  }
  
  ## Load data from .csv file
  gyr_data_raw <- subset(read.csv("data/Gyr_data.csv"), select = -c(1))
  
  ## Rename areas to a more informative name
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
  
  ## Summarise observation by spatial unit and year
  gyr_obs <- gyr_data %>% 
    dplyr::group_by(spatialUnit, Year) %>%
    dplyr::summarise(gyrOccAvg = mean(breeding_attempt, na.rm = TRUE), .groups = "keep")
  
  ## Add year index
  gyr_obs$YearIdx <- gyr_obs$Year - minYear + 1
  
  ## Set up matrix for area-specific data
  gyrOccAvg <- matrix(NA, nrow = N_sUnits, ncol = length(minYear:maxYear))
  
  
  for(x in 1:N_sUnits){
    
    ## Subset data (specific area)
    if(!(sUnits[x] %in% gyr_obs$spatialUnit)){
      stop(paste0("Spatial unit ", sUnits[x], " (index ", x, ") is not in the data."))
    }
    
    gyr_obs_sub <- subset(gyr_obs, spatialUnit == sUnits[x])
    
    ## Extract year-specific data
    for(t in 1:length(minYear:maxYear)){
      
      if(t %in% gyr_obs_sub$YearIdx){
        gyrOccAvg[x, t] <- gyr_obs_sub$gyrOccAvg[which(gyr_obs_sub$YearIdx == t)]
      }
    }
  }
  
  ## Z-standardize covariate values
  meanCov <- mean(gyrOccAvg, na.rm = TRUE)
  sdCov <- sd(gyrOccAvg, na.rm = TRUE)
  gyrOccAvg <- (gyrOccAvg - meanCov) / sdCov
    
  ## Return data
  return(list(gyrOccAvg = gyrOccAvg,
              meanCov = meanCov, 
              sdCov = sdCov))

  }

#d_gyr_occ <- wrangleData_OccGyr(#localities = localities,
#                            areas = areas,
#                            areaAggregation = areaAggregation,
#                            minYear = minYear, maxYear = maxYear)
