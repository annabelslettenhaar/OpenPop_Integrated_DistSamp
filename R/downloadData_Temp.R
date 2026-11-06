#' Download weather data from NVE GridTimeSeries API
#'
#' @param minYear Integer.
#' @param maxYear Integer.
#' @param grid_df Data frame containing grid information.
#' @param param Character string. NVE parameter (default = "tm").
#'
#' @export

downloadData_Temp <- function(minYear,
                                maxYear,
                                grid_df,
                                param = "tm") {
  
  library(httr)
  library(jsonlite)
  library(dplyr)
  
  ## Grid coordinates
  coords_df <- grid_df %>%
    select(grid_id, spatialUnit, x, y) %>%
    mutate(utm_x = round(x), utm_y = round(y), utm_coord = paste(utm_x, utm_y, sep = "/"))
  
  ## Date range
  start.date <- as.Date(paste0(minYear, "-01-01"))
  end.date <- as.Date(paste0(maxYear, "-12-31"))
  
  dates <- seq(start.date, end.date, by = "day")
  
  start_date_str <- format(start.date, "%Y-%m-%d")
  end_date_str <- format(end.date, "%Y-%m-%d")
  
  ## Download
  base_url <- "https://gts.nve.no/api/GridTimeSeries/"
  
  data_list <- list()
  
  for(i in seq_len(nrow(coords_df))){
    coord <- coords_df$utm_coord[i]
    url <- paste0(base_url, coord, "/", start_date_str, "/", end_date_str, "/", param, ".json")
    
    cat("Downloading:", coord, "\n")
    
    res <- httr::GET(url, timeout(seconds = 30))
    dat <- jsonlite::fromJSON(rawToChar(res$content))
    data_list[[coord]] <- dat
  }
  
  df.temp_raw <- bind_rows(
    lapply(names(data_list), function(coord){
      
      data.frame(
        Coordinate = coord,
        Date = dates,
        Value = data_list[[coord]]$Data
      )
    })
  )
  
  return(
    list(
      df.temp_raw = df.temp_raw,
      coords_df = coords_df,
      minYear = minYear,
      maxYear = maxYear
    )
  )
}

temp_raw <- downloadData_Temp(
  minYear = 2013,
  maxYear = 2018,
  grid_df = d_grid$grid_df
)
