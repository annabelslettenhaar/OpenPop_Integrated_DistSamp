#' Extract and prepare rodent abundance covariate data
#'
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

## For testing purposes
#areaAggregation <- TRUE
#areas <- c("Hardangervidda","Dovrefjell", "Børgefjell")
#localities <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#minYear <- 1991
#maxYear <- 2020

wrangleData_RodentGyr <- function(localities = NULL, areas = NULL, areaAggregation, minYear, maxYear){
  
  ## Check if .csv file is available
  if(!file.exists("data/ptar/rodents_1991_2020.csv")){
    stop("Data file (data/ptar/rodents_1991_2020.csv) not found. The workflow currently requires this file as it does not yet support extraction of rodent observation data directly from GBIF/Living Norway.")
  }
  
  ## Load data from .csv file
  rodent_data_raw <- read.csv("data/ptar/rodents_1991_2020.csv")
  
  ## Rename areas
  rodent_data_raw$gyrArea <- as.factor(rodent_data_raw$site)
  rodent_data_raw <- rodent_data_raw %>%
    mutate(gyrArea = recode(gyrArea,
                            "amotsdalen_TOV" = "Hardangervidda",
                            "Mosvatn_TOV" = "Dovrefjell",
                            "Borgefjell_tov" = "Børgefjell"))
  
  ## Filter event data by either locality and year or area and year
  if(areaAggregation){
    rodent_data <- rodent_data_raw %>% 
      dplyr::filter(gyrArea %in% areas) %>%
      dplyr::filter(dplyr::between(year, minYear, maxYear))
  }else{
    rodent_data <- rodent_data_raw %>% 
      dplyr::filter(TerritoryID %in% localities) %>%
      dplyr::filter(dplyr::between(year, minYear, maxYear))
  }
  
  ## Double-check no duplicate records remain
  record_duplicates <- rodent_data %>%
    dplyr::group_by(gyrArea, year) %>%
    dplyr::summarise(recordCount = dplyr::n(), .groups = 'keep') %>%
    dplyr::filter(recordCount > 1)
  
  if(nrow(record_duplicates) > 0){
    stop("There are duplicate transects (> 1 record in same location per year).")
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
    colnames(rodent_data)[which(colnames(rodent_data) == "gyrArea")] <- "spatialUnit"
  }else{
    colnames(rodent_data)[which(colnames(rodent_data) == "locality")] <- "spatialUnit"
  }
  
  ## Create column rodentOcc that balances the observations by trapping effort
  rodent_data <- rodent_data %>%
    mutate(rodentOcc = tot_lem_vol)
  
  ## Summarise observation by spatial unit and year
  rodent_obs <- rodent_data %>% 
    dplyr::group_by(spatialUnit, year) %>%
    dplyr::summarise(rodentAvg = mean(rodentOcc), .groups = "keep")
  
  ## Add year index
  rodent_obs$YearIdx <- rodent_obs$year - minYear + 1
  
  ## Set up matrix for area-specific data
  rodentAvg <- matrix(NA, nrow = N_sUnits, ncol = length(minYear:maxYear))
  
  
  for(x in 1:N_sUnits){
    
    ## Subset data (specific area)
    if(!(sUnits[x] %in% rodent_obs$spatialUnit)){
      stop(paste0("Spatial unit ", sUnits[x], " (index ", x, ") is not in the data."))
    }
    
    rodent_obs_sub <- subset(rodent_obs, spatialUnit == sUnits[x])
    
    ## Extract year-specific data
    for(t in 1:length(minYear:maxYear)){
      
      if(t %in% rodent_obs_sub$YearIdx){
        rodentAvg[x, t] <- rodent_obs_sub$rodentAvg[which(rodent_obs_sub$YearIdx == t)]
      }
    }
  }
  
  ## Z-standardize covariate values
  meanCov <- mean(rodentAvg, na.rm = TRUE)
  sdCov <- sd(rodentAvg, na.rm = TRUE)
  rodentAvg <- (rodentAvg - meanCov) / sdCov
  
  ## Return data
  return(list(rodentAvg = rodentAvg,
              meanCov = meanCov, 
              sdCov = sdCov))
  
}

## Again for testing purposes

#d_rodent <- wrangleData_Rodent(#localities = localities,
#                               areas = areas,
#                               areaAggregation = areaAggregation,
#                               minYear = minYear, maxYear = maxYear)
