
# Plotting posterior densities of covariate slopes

library(ggplot2)
library(dplyr)
library(ggridges)

# Convert to matrix
samps <- as.matrix(IDSM.out)

# Grab only the totDens_raw variables
cov_names <- c("betaR.R", "betaPtar.Occ", "betaPtar.Prod", "betaGyr.S", "betaTemp.R")
samps_sel <- samps[, cov_names, drop = FALSE]

# Prepare data
posterior_df <- samps_sel %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Draw") %>%
  mutate(Parameter = factor(Parameter, levels = rev(unique(Parameter))),
         Parameter = fct_recode(Parameter, 
                            "β-Rodent" = "betaR.R",
                            "β-Occ" = "betaPtar.Occ",
                            "β-Prod" = "betaPtar.Prod",
                            "β-Gyr" = "betaGyr.S",
                            "β-Temp" = "betaTemp.R"))

# Ridge plot with x-axis gradient fill
ggplot(posterior_df, aes(x = Draw, y = Parameter, group = Parameter)) +
  geom_density_ridges_gradient(
    aes(height = after_stat(density), fill = after_stat(x)), # fill by x-axis
    scale = 1,
    rel_min_height = 0.01,
    color = "black"
  ) +
  scale_fill_gradient2(
    low = "#DC0085",
    #mid = "#E4E0DD",
    high = "#11C638",
    midpoint = 0
    #values = c(0, 0.5, 1)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    title = "Posterior distributions of covariate slopes",
    x = "Effect size (per 1 SD increase in predictor)",
    y = " ",
    fill = "Effect size"
  ) +
  theme_minimal(base_size = 14)
