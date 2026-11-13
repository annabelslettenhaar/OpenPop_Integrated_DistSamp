#' Download weather data from NVE GridTimeSeries API
#' 
#' Sample 9 points within each grid cell, extract temperature information from each and then calculate the average
#'
#' @param minYear Integer.
#' @param maxYear Integer.
#' @param grid_df Data frame containing grid information.
#' @param param Character string. NVE parameter (default = "tm").
#'
#' @export

downloadData_Temp <- function(minYear,
                                maxYear,
                                grid_poly,
                                param = "tm") {
  
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(sf)
  
  # ## Grid coordinates
  # coords_df <- grid_df %>%
  #   select(grid_id, spatialUnit, x, y) %>%
  #   mutate(utm_x = round(x), utm_y = round(y), utm_coord = paste(utm_x, utm_y, sep = "/"))
  
  ## Generate 9 points within each grid polygon
  sample_points <- list()
  
  for(i in 1:nrow(grid_poly)){
    poly <- grid_poly[i, ]
    pts <- st_sample(poly, size = 9, type = "regular")
    temp <- st_coordinates(pts)
    sample_points[[i]] <- data.frame(
      grid_id = poly$grid_id,
      spatialUnit = poly$spatialUnit,
      pointID = seq_len(nrow(temp)),
      utm_x = round(temp[,1]),
      utm_y = round(temp[,2]))
  }
  
  coords_df <- bind_rows(sample_points)
  coords_df <- coords_df %>% 
    mutate(utm_coord = paste(utm_x, utm_y, sep="/"))
  
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
    
    # cat("Downloading:", coord, "\n")
    
    res <- httr::GET(url, timeout(seconds = 30))
    dat <- jsonlite::fromJSON(rawToChar(res$content))
    data_list[[coord]] <- dat
  }
  
  coord_lookup <- coords_df %>%
    select(grid_id, spatialUnit, pointID, utm_coord)
  
  df.temp_raw <- bind_rows(
    lapply(names(data_list), function(coord){
      
      idx <- match(coord, coord_lookup$utm_coord)
      
      data.frame(
        grid_id = coord_lookup$grid_id[idx],
        spatialUnit = coord_lookup$spatialUnit[idx],
        pointID = coord_lookup$pointID[idx],
        Coordinate = coord,
        Date = dates,
        Value = data_list[[coord]]$Data
      )
    })
  )
  
  df.temp <- df.temp_raw %>%
    group_by(grid_id, spatialUnit, Date) %>%
    summarise(
      Value = mean(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(
    list(
      df.temp_raw = df.temp,
      coords_df = coords_df,
      minYear = minYear,
      maxYear = maxYear
    )
  )
}

# temp_raw <- downloadData_Temp(
#   minYear = 2013,
#   maxYear = 2018,
#   grid_df = d_grid$grid_df
# )
