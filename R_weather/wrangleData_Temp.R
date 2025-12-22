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
  
wrangleData_Temp <- function(minYear, maxYear, areas, startday, endday) {
  
  library(purrr)
  library(httr)
  library(jsonlite)
  
  # 1. Read ptarmigan transect data
  ptar <- read.csv("data/ptar/event_total.csv")
  
  # 2. Prepare unique transects
  loc <- ptar %>%
    select(footprintWKT, gyrArea)
  
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
  coords_df <- coords_df %>%
    distinct(utm_coord, .keep_all = TRUE)
  
  # 6. Prepare date range
  start.date <- as.Date(paste0(minYear, "-01-01"))
  end.date <- as.Date(paste0(maxYear, "-12-31"))
  dates <- seq(start.date, end.date, "days")
  start_date_str <- format(start.date, "%Y-%m-%d")
  end_date_str <- format(end.date, "%Y-%m-%d")
  
  # 7. Fetch temperature data
  param <- "tm"
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
  
  df.temp <- bind_rows(data_frames)
  
  # 9. Join with area info
  df.temp <- df.temp %>%
    left_join(coords_df[, c("utm_coord", "gyrArea")], by = c("Coordinate" = "utm_coord"))
  
  # 10. Filter chick period (24 June - 15 July)
  df.temp$julianday <- as.numeric(format(df.temp$Date, "%j"))
  springtemp <- subset(df.temp, julianday >= startday & julianday <= endday)
  springtemp$Year <- as.numeric(format(springtemp$Date, "%Y"))
  
  # 11. Summarize mean temperature per area/year
  springtemp_summary <- springtemp %>%
    group_by(gyrArea, Year) %>%
    summarise(temp_mean = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(YearIdx = Year - minYear + 1)
  
  # 12. Create matrix
  sUnits <- areas
  N_sUnits <- length(sUnits)
  mat <- matrix(NA, nrow = N_sUnits, ncol = length(minYear:maxYear))
  
  for (x in seq_len(N_sUnits)) {
    unit <- sUnits[x]
    df_sub <- springtemp_summary[springtemp_summary$gyrArea == unit, ]
    for (t in seq_len(ncol(mat))) {
      if (t %in% df_sub$YearIdx) {
        mat[x, t] <- df_sub$temp_mean[df_sub$YearIdx == t]
      }
    }
  }
  
  # 13. Standardize
  cov_mean <- mean(mat, na.rm = TRUE)
  cov_sd <- sd(mat, na.rm = TRUE)
  mat_std <- (mat - cov_mean) / cov_sd
  
  return(list(
    data = mat_std,
    mean = cov_mean,
    sd = cov_sd
  ))
}
  
# d_temppre <- wrangleData_Temp(minYear = minYear,
#                              maxYear = maxYear,
#                              areas = areas,
#                              startday = 121,
#                              endday = 153)

