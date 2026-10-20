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

wrangleData_Temp <- function(minYear, maxYear, grid_df, startday, endday) {
  
  library(purrr)
  library(httr)
  library(jsonlite)
  
  ## Getting grid coordinates
  coords_df <- grid_df %>%
    select(grid_id, spatialUnit, x, y) %>%
    mutate(utm_x = round(x), utm_y = round(y), utm_coord = paste(utm_x, utm_y, sep = "/"))
  
  ## Set date range (we only want the selected time period)
  start.date <- as.Date(paste0(minYear, "-01-01"))
  end.date <- as.Date(paste0(maxYear, "-12-31"))
  dates <- seq(start.date, end.date, "days")
  start_date_str <- format(start.date, "%Y-%m-%d")
  end_date_str <- format(end.date, "%Y-%m-%d")
  
  ## Use API to download weather data
  param <- "tm"
  base_url <- "http://gts.nve.no/api/GridTimeSeries/"
  coordinates <- coords_df$utm_coord
  data_list <- list()
  
  for(i in seq_len(nrow(coords_df))){
    coord <- coords_df$utm_coord[i]
    url <- paste0(base_url, coord, "/", start_date_str, "/", end_date_str, "/", param, ".json")
    res <- httr::GET(url, timeout(seconds = 30))
    dat <- jsonlite::fromJSON(rawToChar(res$content))
    data_list[[coord]] <- dat
  }
  
  ## Convert to dataframe
  df.temp <- bind_rows(lapply(names(data_list), function(coord){
    data.frame(Coordinate = coord, Date = dates, Value = data_list[[coord]]$Data)
    }))
  
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
  
  ## Setup as matrix for easy input in nimble
  N_grid <- nrow(grid_df)
  N_years <- length(minYear:maxYear)
  
  mat <- matrix(NA, nrow = N_grid, ncol = N_years)
  
  for(g in seq_len(N_grid)){
    this_grid <- grid_df$grid_id[g]
    temp_sub <- springtemp_summary[springtemp_summary$grid_id == this_grid, ]
    for(t in seq_len(N_years)){
      if(t %in% temp_sub$YearIdx){
        mat[g,t] <- temp_sub$temp_mean[temp_sub$YearIdx == t]
      }
    }
  }
  
  ## Standardizing uitput
  cov_mean <- mean(mat, na.rm = TRUE)
  cov_sd <- sd(mat, na.rm = TRUE)
  mat_std <- (mat - cov_mean) / cov_sd
  return(
    list(data = mat_std, mean = cov_mean, sd = cov_sd, grid_lookup = grid_df %>%
        select(grid_id, spatialUnit))
  )
}

# d_temp <- wrangleData_Temp(
#   minYear = minYear,
#   maxYear = maxYear,
#   grid_df = d_grid$grid_df,
#   startday = 121,
#   endday = 153
# )
