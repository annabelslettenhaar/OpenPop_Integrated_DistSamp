#' Extract Weather Data for Ptarmigan Transects
#'
#' This function processes ptarmigan transect data, samples points along transects,
#' retrieves weather data from the NVE GridTimeSeries API for a specified date range
#' and parameter, and calculates mean, standard deviation, and z-standardized values
#' for the chick-rearing period (24 June–15 July) per area and year.
#'
#' @param minYear Integer. The starting year for weather data extraction (e.g., 1990).
#' @param maxYear Integer. The ending year for weather data extraction (e.g., 2020).
#' @param areas string or vector of strings. Names of areas to extract
#' data for.
#' @param startday Integer. Julian day at which the relevant temperature period starts
#' @param endday Integer. Julian day at which the relevant temperature period ends
#'
#' @import purrr httr jsonlite
#' @export

wrangleData_Temp <- function(minYear, maxYear, areas, grid_df, startday, endday) {
  
  # library(purrr)
  # library(httr)
  # library(jsonlite)
  # 
  # ## Getting grid coordinates
  # coords_df <- grid_df %>%
  #   select(grid_id, spatialUnit, x, y) %>%
  #   mutate(utm_x = round(x), utm_y = round(y), utm_coord = paste(utm_x, utm_y, sep = "/"))
  # 
  # ## Set date range (we only want the selected time period)
  # start.date <- as.Date(paste0(minYear, "-01-01"))
  # end.date <- as.Date(paste0(maxYear, "-12-31"))
  # dates <- seq(start.date, end.date, "days")
  # start_date_str <- format(start.date, "%Y-%m-%d")
  # end_date_str <- format(end.date, "%Y-%m-%d")
  # 
  # ## Use API to download weather data
  # param <- "tm"
  # base_url <- "https://gts.nve.no/api/GridTimeSeries/"
  # coordinates <- coords_df$utm_coord
  # data_list <- list()
  # 
  # for(i in seq_len(nrow(coords_df))){
  #   coord <- coords_df$utm_coord[i]
  #   url <- paste0(base_url, coord, "/", start_date_str, "/", end_date_str, "/", param, ".json")
  #   res <- httr::GET(url, timeout(seconds = 30))
  #   dat <- jsonlite::fromJSON(rawToChar(res$content))
  #   data_list[[coord]] <- dat
  # }
  
  data_list <- d_temp_raw$data_list
  coords_df <- d_temp_raw$coords_df
  dates <- d_temp_raw$dates
  minYear <- d_temp_raw$minYear
  maxYear <- d_temp_raw$maxYear
  df.temp <- d_temp_raw$df.temp_raw
  
  ## Join with grid information
  df.temp <- df.temp %>%
    left_join(coords_df %>%
        select(grid_id, spatialUnit, utm_coord), by = c("Coordinate" = "utm_coord"))
  
  ## Select only the spring/early chick period
  df.temp$julianday <- as.numeric(format(df.temp$Date,"%j"))
  springtemp <- df.temp %>%
    filter(julianday >= startday, julianday <= endday)
  springtemp$Year <- as.numeric(format(springtemp$Date,"%Y"))
  
  ## Get average temperature per grid cell and year
  springtemp_summary <- springtemp %>%
    group_by(grid_id, Year) %>%
    summarise(temp_mean = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(YearIdx = Year - minYear + 1)
  
  ## Setup array x, g, t for easy input in nimble
  areas <- areas
  N_areas <- length(areas)
  grid_count <- grid_df %>%
    dplyr::group_by(spatialUnit) %>%
    dplyr::summarise(count = n())
  N_grid <- grid_count$count
  N_years <- length(minYear:maxYear)
  temp_array <- array(0, dim = c(N_areas, max(N_grid), N_years))
  
  for(x in 1:N_areas){
    grid_area <- grid_df[grid_df$spatialUnit == areas[x], ]
    grid_area <- dplyr::arrange(grid_area, grid_id)
    N_grid_x <- nrow(grid_area)
    for(g in 1:N_grid_x){
      this_grid <- grid_area$grid_id[g]
      temp_sub <- springtemp_summary[springtemp_summary$grid_id == this_grid, ]
      for(t in 1:N_years){
        if(t %in% temp_sub$YearIdx){
          temp_array[x,g,t] <- temp_sub$temp_mean[temp_sub$YearIdx == t]
        }
      }
    }
  }
  
  ## Long-term mean temperature per grid cell
  mean_temp_array <- apply(temp_array, c(1,2), mean, na.rm = TRUE)
  
  ## Temperature anomalies
  temp_anom_array <- array(NA, dim = dim(temp_array))
  
  for(x in 1:N_areas){
    for(g in 1:max(N_grid)){
      temp_anom_array[x,g,] <- temp_array[x,g,] - mean_temp_array[x,g]
    }
  }
  
  ## Standardizing uitput (separately for mean and anomaly)
  mean_temp_mean <- mean(mean_temp_array, na.rm = TRUE)
  mean_temp_sd <- sd(mean_temp_array, na.rm = TRUE)
  mean_temp_std <- (mean_temp_array - mean_temp_mean) / mean_temp_sd
  
  anom_mean <- mean(temp_anom_array, na.rm = TRUE)
  anom_sd <- sd(temp_anom_array, na.rm = TRUE)
  temp_anom_std <- (temp_anom_array - anom_mean) / anom_sd
  
  return(
    list(
      MeanTemp = mean_temp_std,
      TempAnom = temp_anom_std,
      meanTemp_mean = mean_temp_mean,
      meanTemp_sd = mean_temp_sd,
      anom_mean = anom_mean,
      anom_sd = anom_sd,
      N_grid = N_grid,
      areas = areas))
}

# d_temp <- wrangleData_Temp(
#   minYear = minYear,
#   maxYear = maxYear,
#   grid_df = d_grid$grid_df,
#   startday = 121,
#   endday = 153
# )
