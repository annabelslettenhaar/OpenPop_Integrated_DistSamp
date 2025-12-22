library(coda)
library(dplyr)
library(stringr)
library(tidyr)

posterior <- as.matrix(temprun)

# keep only the parameters of interest
effects <- posterior[, grep("betaPtar.Prod|betaPtar.Occ|betaGyr.S|betaR.R|betaTemp.R", colnames(posterior))]

# reshape to long format
effects_df <- as.data.frame(effects) %>%
  tibble::rownames_to_column("iteration") %>%
  pivot_longer(-iteration, names_to = "parameter", values_to = "value") %>%
  mutate(
    type = case_when(
      str_detect(parameter, "betaPtar.Occ") ~ "Prey → Predator occupancy",
      str_detect(parameter, "betaPtar.Prod") ~ "Prey → Predator productivity",
      str_detect(parameter, "betaGyr.S") ~ "Predator → Prey surival",
      str_detect(parameter, "betaR.R") ~ "Alternative Prey → Prey recruitment",
      str_detect(parameter, "betaTemp.R") ~ "Temperature → Prey recruitment"
    ),
    area = as.integer(str_extract(parameter, "(?<=\\[)\\d+(?=\\])"))
  )

summary_df <- effects_df %>%
  group_by(type, area) %>%
  summarise(
    median = median(value),
    lci = quantile(value, 0.025),
    uci = quantile(value, 0.975),
    .groups = "drop"
  )

area_names <- c(
  "1" = "Hardangervidda",
  "2" = "Dovrefjell",
  "3" = "Børgefjell")
summary_df <- summary_df %>%
  mutate(area_name = recode(as.character(area), !!!area_names))

ggplot(summary_df, aes(x = median, y = factor(area_name), color = type)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = lci, xmax = uci),
                 height = 0.2,
                 position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  theme(text = element_text(size=16)) +
  labs(x = "Posterior effect (median ± 95% CI)", y = "Area", color = "Effect")


