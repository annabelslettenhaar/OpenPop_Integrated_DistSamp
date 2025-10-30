#' Extract Weather Data for Ptarmigan Transects
#'
#' This function processes ptarmigan transect data, samples points along transects,
#' retrieves snow depth data from the NVE GridTimeSeries API for a specified date 
#' range and parameter, and calculates for every year at which julian day the first 
#' 7-day snow-free period starts. 
#'
#' @param minYear Integer. The starting year for weather data extraction (e.g., 1990).
#' @param maxYear Integer. The ending year for weather data extraction (e.g., 2020).
#' @param areas string or vector of strings. Names of areas to extract
#' data for.
#'
#' @import purrr httr jsonlite sf zoo
#' @export

wrangleData_Snow <- function(minYear, maxYear, areas) {
  library(purrr)
  library(httr)
  library(jsonlite)
  library(sf)
  library(zoo)
  
  # 1. Read ptarmigan transect data
  ptar <- read.csv("data/ptar/event_total.csv")
  
  # 2. Prepare unique transects
  loc <- ptar %>% select(footprintWKT, gyrArea)
  unique_coord <- loc %>%
    distinct(footprintWKT, .keep_all = TRUE) %>%
    filter(!is.na(footprintWKT))
  
  # 3. Convert WKT to sf and transform to UTM
  transects <- st_as_sf(unique_coord, wkt = "footprintWKT", crs = 4326)
  transects_utm <- st_transform(transects, 32633)
  
  # 4. Sample points every 500 m
  sampled_points <- map_dfr(1:nrow(transects_utm), function(i) {
    line <- transects_utm[i, ]
    len <- as.numeric(st_length(line))
    step <- 500 / len
    sample_seq <- unique(c(seq(0, 1, by = step), 1))
    pts <- st_line_sample(line$footprintWKT, sample = sample_seq)
    pts_sf <- st_sf(
      gyrArea = line$gyrArea,
      geometry = st_cast(pts, "POINT"),
      crs = st_crs(line)
    )
    return(pts_sf)
  })
  
  # 5. Extract UTM coordinates
  utm_coords <- as.data.frame(st_coordinates(sampled_points))
  coords_df <- cbind(st_drop_geometry(sampled_points), utm_coords)
  coords_df$utm_x <- round(coords_df$X)
  coords_df$utm_y <- round(coords_df$Y)
  coords_df$utm_coord <- as.factor(paste(coords_df$utm_x, coords_df$utm_y, sep="/"))
  coords_df <- coords_df %>% distinct(utm_coord, .keep_all = TRUE)
  
  # 6. Prepare date range
  start.date <- as.Date(paste0(minYear, "-01-01"))
  end.date <- as.Date(paste0(maxYear, "-12-31"))
  dates <- seq(start.date, end.date, "days")
  start_date_str <- format(start.date, "%Y-%m-%d")
  end_date_str <- format(end.date, "%Y-%m-%d")
  
  # 7. Fetch snow depth data
  param <- "sd"
  base_url <- "http://gts.nve.no/api/GridTimeSeries/"
  coordinates <- coords_df$utm_coord
  data_list <- list()
  
  for (coord in coordinates) {
    url <- paste0(base_url, coord, "/", start_date_str, "/", end_date_str, "/", param, ".json")
    res <- httr::GET(url, timeout(seconds = 30))
    dat <- jsonlite::fromJSON(rawToChar(res$content))
    data_list[[coord]] <- dat
  }
  
  # 8. Combine into dataframe
  data_frames <- list()
  for (coord_key in names(data_list)) {
    dat <- data_list[[coord_key]]
    data_frame <- data.frame(Coordinate = coord_key,
                             Date = dates,
                             Value = dat$Data)
    data_frames[[coord_key]] <- data_frame
  }
  
  df.snow <- bind_rows(data_frames)
  
  # 9. Join with area info
  df.snow <- df.snow %>%
    left_join(coords_df[, c("utm_coord", "gyrArea")], by = c("Coordinate" = "utm_coord"))
  
  # 10. Add Julian day and year
  df.snow$julianday <- as.numeric(format(df.snow$Date, "%j"))
  df.snow$Year <- as.numeric(format(df.snow$Date, "%Y"))
  
  # 11. Detect first snow-free week
  snowfree_week <- df.snow %>%
    arrange(gyrArea, Year, julianday) %>%
    group_by(gyrArea, Year) %>%
    mutate(snow_free = Value == 0) %>%
    mutate(consecutive_free_week = zoo::rollapply(snow_free, 7, all, fill = NA, align = "left")) %>%
    filter(consecutive_free_week == TRUE) %>%
    summarise(snowfree_week_start = min(julianday), .groups = "drop") %>%
    mutate(YearIdx = Year - minYear + 1)
  
  # 12. Create matrix
  sUnits <- areas
  N_sUnits <- length(sUnits)
  mat <- matrix(NA, nrow = N_sUnits, ncol = length(minYear:maxYear))
  
  for (x in seq_len(N_sUnits)) {
    unit <- sUnits[x]
    df_sub <- snowfree_week[snowfree_week$gyrArea == unit, ]
    for (t in seq_len(ncol(mat))) {
      if (t %in% df_sub$YearIdx) {
        mat[x, t] <- df_sub$snowfree_week_start[df_sub$YearIdx == t]
      }
    }
  }
  
  # 14. Standardize the matrix
  cov_mean <- mean(mat, na.rm = TRUE)
  cov_sd <- sd(mat, na.rm = TRUE)
  mat_std <- (mat - cov_mean) / cov_sd
  
  # 15. Return both raw and standardized data
  return(list(
    #average_snow_depth_may20 = mat,
    standardized = mat_std,
    mean = cov_mean,
    sd = cov_sd
  ))
  return(list(
    data = mat
  ))
}

# d_snow <- wrangleData_SnowFreeWeek(minYear = minYear, 
#                                    maxYear = maxYear,
#                                    areas = areas)
