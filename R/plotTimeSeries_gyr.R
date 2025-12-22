#' Plot Time Series of Posterior Summaries
#'
#' Creates two ggplot visualizations from posterior samples:  
#' (1) Faceted plots for ptarmigan density, gyrfalcon pressure, and productivity.  
#' (2) Combined plot of ptarmigan density and gyrfalcon pressure with dual y-axes.
#'
#' @param IDSM_out A matrix or object convertible to a matrix containing posterior samples.
#' @param minYear An integer specifying the first year in the time series.
#'
#' @return A list with two `ggplot` objects:  
#' \itemize{
#'   \item `p_prey`: Faceted plot of all three variables.
#'   \item `p_prey_comb`: Combined plot of ptarmigan density and gyrfalcon pressure.
#' }
#' @import ggplot2 cowplot dplyr tidyr stringr tibble
#' @export

plotTimeSeries_gyr <- function(IDSM_out, minYear) {
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  
  # Convert to matrix and long format
  samps <- as.matrix(IDSM_out)
  long_df <- as.data.frame(samps) %>%
    tibble::rownames_to_column("iteration") %>%
    pivot_longer(-iteration, names_to = "param", values_to = "value") %>%
    filter(str_detect(param, "totDens_raw|GyrPressure_raw|terrProd"))
  
  # Extract variable type, area, and year
  parsed <- str_match(long_df$param, "^([a-zA-Z_]+)\\[(\\d+),\\s*(\\d+)\\]$")
  long_df$variable <- parsed[, 2]
  long_df$area_index <- as.integer(parsed[, 3])
  long_df$year_index <- as.integer(parsed[, 4])
  
  # Map year/area index to actual values
  long_df <- long_df %>%
    mutate(year = minYear + year_index - 1,
           area = recode(area_index,
                         `1` = "Hardangervidda",
                         `2` = "Dovrefjell",
                         `3` = "Børgefjell"),
           value = if_else(variable == "totDens_raw", value * 1e6, value)) %>%
    filter(variable %in% c("totDens_raw","GyrPressure_raw","terrProd"))
  
  # Summarize posterior samples
  summary_df <- long_df %>%
    group_by(variable, area, year) %>%
    summarise(
      mean = mean(value),
      lower = quantile(value, 0.025),
      upper = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(Variable = recode(variable,
                             "GyrPressure_raw" = "GyrPressure",
                             "terrProd" = "Gyrfalcon productivity",
                             "totDens_raw" = "Ptarmigan density (birds/km^2)"))
  
  summary_df$area <- factor(summary_df$area, levels = c("Hardangervidda", "Dovrefjell", "Børgefjell"))
  
  ## Plot 1: All variables faceted
  p_prey <- ggplot(summary_df, aes(x = year, y = mean, color = Variable)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Variable),
                alpha = 0.2, color = NA) +
    facet_grid(Variable ~ area, scales = "free_y") +
    scale_color_manual(values = c(
      "Ptarmigan density (birds/km^2)" = "#7EBB4C",
      "GyrPressure" = "#E05D00",
      "Gyrfalcon productivity" = "#CD348E"
    )) +
    scale_fill_manual(values = c(
      "Ptarmigan density (birds/km^2)" = "#7EBB4C",
      "GyrPressure" = "#E05D00",
      "Gyrfalcon productivity" = "#CD348E"
    )) +
    theme_minimal(base_size = 14) +
    theme(strip.text.y = element_blank(),
          legend.position = "bottom")
  
  ## Plot 2: Combined ptarmigan density & gyrpressure
  scale_factor <- max(summary_df$mean[summary_df$Variable == "Ptarmigan density (birds/km^2)"]) /
    max(summary_df$mean[summary_df$Variable == "GyrPressure"])
  
  p_prey_comb <- ggplot(summary_df, aes(x = year)) +
    geom_line(data = subset(summary_df, Variable == "Ptarmigan density (birds/km^2)"),
              aes(y = mean, color = "Ptarmigan density (birds/km^2)")) +
    geom_ribbon(data = subset(summary_df, Variable == "Ptarmigan density (birds/km^2)"),
                aes(ymin = lower, ymax = upper, fill = "Ptarmigan density (birds/km^2)"),
                alpha = 0.2, color = NA) +
    geom_line(data = subset(summary_df, Variable == "GyrPressure"),
              aes(y = mean * scale_factor, color = "GyrPressure")) +
    geom_ribbon(data = subset(summary_df, Variable == "GyrPressure"),
                aes(ymin = lower * scale_factor, ymax = upper * scale_factor, fill = "GyrPressure"),
                alpha = 0.2, color = NA) +
    facet_wrap(~ area, scales = "free_y") +
    scale_y_continuous(
      name = "Ptarmigan density (birds/km²)",
      sec.axis = sec_axis(~ . / scale_factor, name = "GyrPressure")
    ) +
    scale_color_manual(values = c(
      "Ptarmigan density (birds/km^2)" = "#7EBB4C",
      "GyrPressure" = "#E05D00"
    )) +
    scale_fill_manual(values = c(
      "Ptarmigan density (birds/km^2)" = "#7EBB4C",
      "GyrPressure" = "#E05D00"
    )) +
    guides(fill = "none") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
                   strip.text = element_text(face = "bold")) +
    labs(color = "Variable", fill = "Variable")
  
  return(list(p_prey = p_prey, p_prey_comb = p_prey_comb))
}

# plotTimeSeries(minYear = minYear,
#                           IDSM_out = IDSM.out)
#   