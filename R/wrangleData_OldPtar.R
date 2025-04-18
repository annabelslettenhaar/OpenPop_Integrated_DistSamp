#' Extract and filter transect and observation data from DwC archive
#'
#' @param localities string or vector of strings. Names of localities to extract
#' data for. Either localities or areas must be provided. 
#' @param areas string or vector of strings. Names of areas to extract
#' data for. Either localities or areas must be provided.
#' @param areaAggregation logical. If TRUE, areas are used as smallest spatial unit. If FALSE, locations (within areas) are used as smallest spatial unit.
#' @param minYear integer. Earliest year of data to extract.
#' @param maxYear integer. Latest year of data to extract.  
#'
#' @return list of 3 tibbles. `d_trans` contains information on transects 
#' (events). `d_obs`contains information on observations made along transects 
#' (distance to transect line, numbers of birds in each age/sex class observed,
#' etc.). `d_coord` contains averaged coordinates for the selected spatial units. 
#' 
#' @export
#'
#' @examples

## For testing purposes
#areaAggregation <- TRUE
#areas <- c("Åmotsdalen", "Børgefjell", "Møsvatn")
#localities <- c("Gåvålia", "TOV-Åmotsdalen", "TOV-Børgefjell", "TOV-Møsvatn")
#minYear <- 1991
#maxYear <- 2014


wrangleData_OldPtar <- function(localities = NULL, areas = NULL, areaAggregation, minYear, maxYear){
  
  ## Check if .csv files are available
  files <- c("data/ptar/event_1991_2014.csv", "data/ptar/occurrence_1991_2014.csv")
  
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop(paste("The following required files are missing:", paste(missing_files, collapse = ", ")))
  }
  
  ## Load data from .csv files
  Eve <- read.csv("data/ptar/event_1991_2014.csv")
  Occ <- read.csv("data/ptar/occurrence_1991_2014.csv")
  
  ## Filter event data by either locality and year or area and year
  if(areaAggregation){
    Eve <- Eve %>% 
      dplyr::filter(area %in% areas) %>%
      dplyr::filter(dplyr::between(year, minYear, maxYear))
  }else{
    Eve <- Eve %>% 
      dplyr::filter(dplyr::between(year, minYear, maxYear))
  }
  
  ## Assemble transect level info
  d_trans <- Eve %>% 
    dplyr::select(locationID, eventDate, eventID, eventRemarks, transectlength, 
                  stateProvince, municipality, area, year) %>%
    dplyr::filter(eventRemarks == "Line transect")
  
  ## Identify and remove transects with suspiciously short length (< 200 m) and duplicate transects
  bad_transects <- d_trans$eventID[which(d_trans$transectlength < 200)]
  d_trans <- d_trans %>%
    dplyr::filter(!(eventID %in% bad_transects))
  
  ## Double-check no duplicate transects remain
  transect_duplicates <- d_trans %>%
    dplyr::group_by(locationID, area, year) %>%
    dplyr::summarise(transectCount = dplyr::n(), .groups = 'keep') %>%
    dplyr::filter(transectCount > 1)
  
  if(nrow(transect_duplicates) > 0){
    stop("There are duplicate transects (> 1 transect in same location per year).")
  }
  
  ## Assemble observation level info
  
  # Observations: distance to transect lines
  d_obsTemp <- Eve %>% 
    dplyr::select(locationID, area, parentEventID, eventID, eventRemarks, 
                  eventDate, year, linedist) %>%
    dplyr::filter(eventRemarks == "Human observation") %>%
    dplyr::rename(DistanceToTransectLine = linedist) %>%
    dplyr::select(locationID, area, parentEventID, eventID, DistanceToTransectLine, year)
  
  # Observations: remaining information (willow ptarmigan only)
  d_obs <- Occ %>% 
    dplyr::filter(!is.na(eventID)) %>%
    dplyr::select(eventID, scientificName, individualCount, sex, lifeStage) %>%
    dplyr::mutate(lifeStage = if_else(is.na(lifeStage), "unknown", lifeStage),
                  sex = if_else(is.na(sex), "unknown", sex)) %>%
    dplyr::mutate(SexStage = str_c(sex, lifeStage, sep = "")) %>%
    dplyr::select(-sex, -lifeStage) %>%
    tidyr::spread(key = "SexStage", value = "individualCount", fill = 0) %>%
    dplyr::right_join(., d_obsTemp, by = "eventID") %>%
    dplyr::filter(scientificName == "Lagopus lagopus")
  
  ## Remove observations from transects with suspiciously short length (< 200) and duplicate transects
  d_obs <- d_obs %>%
    dplyr::filter(!(parentEventID %in% bad_transects))
  
  ## Extract locality-/area-level coodinate information
  
  # Rename appropriate column
  if(areaAggregation){
    colnames(Eve)[which(colnames(Eve) == "area")] <- "spatialUnit"
  }else{
    colnames(Eve)[which(colnames(Eve) == "locality")] <- "spatialUnit"
  }
  
  # Calculate average longitude and latitude per locality/area
  d_coord <- Eve %>%
    dplyr::select("locationID", "spatialUnit") %>%
    dplyr::distinct() %>%
    dplyr::group_by(spatialUnit)
  
  ## Collate and return data
  LT_data <- list(d_trans = d_trans, d_obs = d_obs, d_coord = d_coord)
}

#LT_data <- wrangleData_DwCPtar(#localities = localities,
#                               areas = areas,
#                               areaAggregation = areaAggregation,
#                               minYear = minYear, maxYear = maxYear)
