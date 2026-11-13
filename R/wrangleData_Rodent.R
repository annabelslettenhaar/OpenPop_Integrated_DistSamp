#' Extract and prepare rodent abundance covariate data
#'
#' @param duplTransects vector of strings. (parent) EventIDs of transects that are duplicates and need removing from the data
#' @param localities string or vector of strings. Names of localities to extract
#' data for. Either localities or areas must be provided. 
#' @param areas string or vector of strings. Names of areas to extract
#' data for. Either localities or areas must be provided.
#' @param areaAggregation logical. If TRUE, areas are used as smallest spatial unit. If FALSE, locations (within areas) are used as smallest spatial unit.
#' @param minYear integer. Earliest year of data to extract.
#' @param maxYear integer. Latest year of data to extract.  
#'
#' @return a matrix containing the average number of transects with rodent observations per area and year.
#' @export
#'
#' @examples

wrangleData_Rodent <- function(duplTransects, localities = NULL, areas = NULL, areaAggregation, minYear, maxYear){
  
  ## Check if .rds file is available
  if(!file.exists("data/Rodent_data.rds")){
    stop("Data file (data/Rodent_data.rds) not found. The workflow currently requires this file as it does not yet support extraction of rodent observation data directly from GBIF/Living Norway.")
  }
  
  ## Load data from .rds file
  rodent_data_raw <- readRDS("data/Rodent_data.rds")
  
  ## Filter event data by either locality and year or area and year
  if(areaAggregation){
    rodent_data <- rodent_data_raw %>% 
      dplyr::mutate(date = as.Date(date)) %>%
      dplyr::mutate(Year = lubridate::year(date)) %>%
      dplyr::filter(verbatimLocality %in% areas) %>%
      dplyr::filter(dplyr::between(Year, minYear, maxYear))
  }else{
    rodent_data <- rodent_data_raw %>% 
      dplyr::mutate(date = as.Date(date)) %>%
      dplyr::mutate(Year = lubridate::year(date)) %>%
      dplyr::filter(locality %in% localities) %>%
      dplyr::filter(dplyr::between(Year, minYear, maxYear))
  }
  
  ## Identify and remove transects with suspiciously short length (< 200 m) and duplicate transects
  bad_transects <- c(rodent_data$EventID[which(rodent_data$sampleSizeValue < 200)], duplTransects)
  rodent_data <- rodent_data %>%
    dplyr::filter(!(EventID %in% bad_transects))
  
  ## Double-check no duplicate transects remain
  transect_duplicates <- rodent_data %>%
    dplyr::group_by(locationID, locality, verbatimLocality, Year) %>%
    dplyr::summarise(transectCount = dplyr::n(), .groups = 'keep') %>%
    dplyr::filter(transectCount > 1)
  
  if(nrow(transect_duplicates) > 0){
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
    colnames(rodent_data)[which(colnames(rodent_data) == "verbatimLocality")] <- "spatialUnit"
  }else{
    colnames(rodent_data)[which(colnames(rodent_data) == "locality")] <- "spatialUnit"
  }
  
  ## Merge with grid cell information/lookup table
  rodent_data <- rodent_data %>%
    left_join(
      d_grid$grid_area_lookup %>%
        select(locationID, grid_id),
      by = "locationID",
      relationship = "many-to-many"
    )
  
  if(any(is.na(rodent_data$grid_id))){
    stop("Some locationID values could not be matched to a grid_id.")
  }
  
  # ## Summarise observation by spatial unit and year
  # rodent_obs <- rodent_data %>% 
  #   dplyr::group_by(spatialUnit, Year) %>%
  #   dplyr::summarise(rodentAvg = mean(rodentOcc, na.rm = TRUE),
  #                    .groups = "keep")
  
  
  ## Summarise observation by grid cell, spatial unit, and year
  rodent_obs <- rodent_data %>% 
    dplyr::group_by(grid_id, spatialUnit, Year) %>%
    dplyr::summarise(rodentAvg = mean(rodentOcc, na.rm = TRUE),
                     .groups = "drop")
  
  
  ## Add year index
  rodent_obs$YearIdx <- rodent_obs$Year - minYear + 1
    
  # ## Set up matrix for area-specific data
  # rodentAvg <- matrix(NA, nrow = N_sUnits, ncol = length(minYear:maxYear))

  ## Set up array for gridcell level data
  
  grid_count <- d_grid$grid_df %>%
    dplyr::group_by(spatialUnit) %>%
    dplyr::summarise(count = n())
  
  N_grid <- grid_count$count
  N_years <- length(minYear:maxYear)
  N_areas <- length(areas)
  
  rodent_array <- array(NA, dim = c(N_areas, max(N_grid), N_years))
  
  # for(x in 1:N_sUnits){
  #   
  #   ## Subset data (specific area)
  #   if(!(sUnits[x] %in% rodent_obs$spatialUnit)){
  #     stop(paste0("Spatial unit ", sUnits[x], " (index ", x, ") is not in the data."))
  #   }
  #   
  #   rodent_obs_sub <- subset(rodent_obs, spatialUnit == sUnits[x])
  # 
  #   ## Extract year-specific data
  #   for(t in 1:length(minYear:maxYear)){
  #     
  #     if(t %in% rodent_obs_sub$YearIdx){
  #       rodentAvg[x, t] <- rodent_obs_sub$rodentAvg[which(rodent_obs_sub$YearIdx == t)]
  #     }
  #   }
  # }
  
  ## Fill the array
  
  for(x in 1:N_areas){
    
    grid_area <- d_grid$grid_df[d_grid$grid_df$spatialUnit == areas[x], ]
    grid_area <- dplyr::arrange(grid_area, grid_id)
    N_grid_x <- nrow(grid_area)
    
    for(g in 1:N_grid_x){
      
      this_grid <- grid_area$grid_id[g]
      rodent_sub <- rodent_obs[rodent_obs$grid_id == this_grid & rodent_obs$spatialUnit == areas[x], ]
      
      for(t in 1:N_years){
        
        if(t %in% rodent_sub$YearIdx){
          
          rodent_array[x,g,t] <- rodent_sub$rodentAvg[rodent_sub$YearIdx == t]
        }
      }
    }
  }
  
 
  ## Z-standardize covariate values
  # meanCov <- mean(rodentAvg, na.rm = TRUE)
  # sdCov <- sd(rodentAvg, na.rm = TRUE)
  # rodentAvg <- (rodentAvg - meanCov) / sdCov
  meanCov <- mean(rodent_array, na.rm = TRUE) ## Maybe eventually we want to correct for transect effort (length)
  sdCov <- sd(rodent_array, na.rm = TRUE)
  rodentAvg <- (rodent_array - meanCov) / sdCov
  
  ## Return data
  return(list(rodentAvg = rodentAvg,
              meanCov = meanCov, 
              sdCov = sdCov,
              N_grid = N_grid,
              areas = areas))
  
}