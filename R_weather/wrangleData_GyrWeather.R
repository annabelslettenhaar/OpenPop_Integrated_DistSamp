#' Extract Weather Data from gyrfalcon territory 'cores'
#'
#' The raw data input is extracted from gridded snow and weather data provided by 
#' SeNorge. The extraction process is done locally due to sensitive nest coordinates
#' needed for extraction of these data. 
#' 
#' We calculate mean temperature during territory establishment (February)
#' We calculate mean snow depth during territory establishment
#' We calculate mean temperature during the early breeding phase (egg laying, April)
#' We calculate mean temperature during the chick phase (between 15th of May (135) 
#' to 1st of July (182))
#' We calculate 
#' We calculate mean snow depth during chick development (20th of May)
#' 
#' The output is one z-standardized mean value per area per year. 
#'
#' @param minYear Integer. The starting year for weather data extraction (e.g., 1990).
#' @param maxYear Integer. The ending year for weather data extraction (e.g., 2020).
#' @param areas string or vector of strings. Names of areas to extract
#' data for.
#' @param byArea TRUE/FALSE. TRUE: calculate means per area per year, FALSE: calculate
#' means per year only. 
#'
#' @import 
#' @export

wrangleData_GyrWeather <- function(minYear, maxYear, areas, byArea) {
  
  # Load and preprocess
  raw_weather <- read.csv("data/weather/gyrterr_snow_temp_alldata.csv")
  raw_weather$Date <- as.Date(raw_weather$Date)
  raw_weather$julianday <- as.numeric(format(raw_weather$Date, "%j"))
  raw_weather$year <- as.numeric(format(raw_weather$Date, "%Y"))
  raw_weather$dateonly <- format(raw_weather$Date, "%m-%d")
  
  raw_weather$gyrArea <- as.factor(raw_weather$Area)
  raw_weather <- raw_weather %>%
    mutate(gyrArea = recode(gyrArea,
                            "1" = "Hardangervidda",
                            "2" = "Dovrefjell",
                            "3" = "Børgefjell")) %>%
    filter(year >= minYear & year <= maxYear,
           gyrArea %in% areas)
  
  years <- minYear:maxYear
  
  # Helper to create standardized matrix
  make_matrix <- function(df, value_col) {
    if (byArea) {
      mat <- matrix(NA, nrow = length(areas), ncol = length(years))
      for (i in seq_along(areas)) {
        for (j in seq_along(years)) {
          val <- df[[value_col]][df$gyrArea == areas[i] & df$year == years[j]]
          if (length(val) > 0) {
            mat[i, j] <- val
          }
        }
      }
    } else {
      mat <- matrix(NA, nrow = 1, ncol = length(years))
      for (j in seq_along(years)) {
        val <- df[[value_col]][df$year == years[j]]
        if (length(val) > 0) {
          mat[1, j] <- mean(val, na.rm = TRUE)
        }
      }
    }
    
    cov_mean <- mean(mat, na.rm = TRUE)
    cov_sd <- sd(mat, na.rm = TRUE)
    mat_std <- unname((mat - cov_mean) / cov_sd)
    return(list(data = mat_std, mean = cov_mean, sd = cov_sd))
  }
  
  # Conditional grouping
  group_vars <- if (byArea) c("gyrArea", "year") else "year"
  
  # Chick period
  chick_df <- raw_weather %>%
    filter(julianday >= 135 & julianday <= 182) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(temp_mean = mean(temp, na.rm = TRUE), .groups = "drop")
  
  # April
  april_df <- raw_weather %>%
    filter(month(Date) == 4) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(temp_mean = mean(temp, na.rm = TRUE), .groups = "drop")
  
  # February
  feb_df <- raw_weather %>%
    filter(month(Date) == 2) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(temp_mean = mean(temp, na.rm = TRUE),
              snow_mean = mean(snow, na.rm = TRUE), .groups = "drop")
  
  # Snow depth on 20 May
  sd20_df <- raw_weather %>%
    filter(dateonly == "05-20") %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(snow_mean = mean(snow, na.rm = TRUE), .groups = "drop")
  
  # Total precipitation during nestling period
  precip <- raw_weather %>%
    filter(julianday >= 135 & julianday <= 182) %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(precip_total = sum(precip, na.rm = TRUE), .groups = "drop")
  
  # Number of periods with 5 consecutive rainy days during nestling period
  rain_streaks <- raw_weather %>%
    filter(julianday >= 135 & julianday <= 182) %>%
    group_by(across(all_of(group_vars))) %>%
    arrange(julianday) %>%
    mutate(rainy = precip > 1) %>%
    summarise(
      long_streaks = {
        rle_result <- rle(rainy)
        sum(rle_result$values & rle_result$lengths >= 5)
      },
      .groups = "drop"
    )
  
  
  # Return list of standardized matrices
  return(list(
    chicktemp = make_matrix(chick_df, "temp_mean"),
    apriltemp = make_matrix(april_df, "temp_mean"),
    febtemp = make_matrix(feb_df, "temp_mean"),
    febsnow = make_matrix(feb_df, "snow_mean"),
    SD20 = make_matrix(sd20_df, "snow_mean"),
    precip = make_matrix(precip, "precip_total"),
    longrain = make_matrix(rain_streaks, "long_streaks")
  ))
}


# test <- wrangleData_GyrWeather(minYear = minYear,
#                                maxYear = maxYear,
#                                areas = areas,
#                                byArea = TRUE)
