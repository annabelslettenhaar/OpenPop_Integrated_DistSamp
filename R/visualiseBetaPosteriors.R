#' Plot Posterior Densities of Covariate Slopes
#'
#' Creates a ridge plot of posterior distributions for selected covariate slopes
#' from the IPM model output.
#'
#' @param IDSM_out A matrix or object convertible to a matrix containing posterior samples.
#' @param cov_names A character vector of parameter names to include in the plot.
#'
#' @return A `ggplot` object showing ridge plots of posterior densities.
#'
#' @import ggplot2 dplyr ggridges forcats tidyr
#' @export

visualiseBetaPosteriors <- function(IDSM_out, cov_names) {
  library(ggplot2)
  library(dplyr)
  library(ggridges)
  library(forcats)
  library(tidyr)
  
  # Convert to matrix and select covariates
  samps <- as.matrix(IDSM_out)
  samps_sel <- samps[, cov_names, drop = FALSE]
  
  # Prepare data
  posterior_df <- samps_sel %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Draw") %>%
    mutate(Parameter = factor(Parameter, levels = rev(unique(Parameter))),
           Parameter = fct_recode(Parameter,
                                  "β-Rodent" = "betaR.R",
                                  "β-Init" = "betaPtar.Occ",
                                  "β-Prod" = "betaPtar.Prod",
                                  "β-Gyr" = "betaGyr.S",
                                  "β-Temp" = "betaTemp.R"))
  
  # Ridge plot
  p <- ggplot(posterior_df, aes(x = Draw, y = Parameter, group = Parameter)) +
    geom_density_ridges_gradient(
      aes(height = after_stat(density), fill = after_stat(x)),
      scale = 1,
      rel_min_height = 0.01,
      color = "black"
    ) +
    scale_fill_gradient2(
      low = "#DC0085",
      high = "#11C638",
      midpoint = 0
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = "Posterior distributions of covariate slopes",
      x = "Effect size (per 1 SD increase in predictor)",
      y = NULL,
      fill = "Effect size"
    ) +
    theme_minimal(base_size = 14)
  
  return(p)
}

# visualiseBetaPosteriors(IDSM_out = IDSM.out,
#                           cov_names = c("betaR.R", "betaPtar.Occ",
#                                         "betaPtar.Prod", "betaGyr.S",
#                                         "betaTemp.R"))
# 

