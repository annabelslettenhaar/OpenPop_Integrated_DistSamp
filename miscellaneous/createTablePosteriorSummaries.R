
## Make a table with the posterior summaries

library(flextable)

survival_bl <- as.data.frame(postSum_list$Mu.S)
survival_bl <- survival_bl %>%
  select(area, Median, uCI, lCI, SD, CV) %>%
  mutate(Species = "Ptarmigan",
         Variable = "μ_S",
         Definition = "Baseline survival probability")

survival <- as.data.frame(postSum_list$S)
survival <- survival %>%
  select(area, Median, uCI, lCI, SD, CV) %>%
  mutate(Species = "Ptarmigan",
         Variable = "S",
         Definition = "Survival probability")

recruitment_bl <- as.data.frame(postSum_list$Mu.R)
recruitment_bl <- recruitment_bl %>%
  select(area, Median, uCI, lCI, SD, CV) %>%
  mutate(Species = "Ptarmigan",
         Variable = "μ_R",
         Definition = "Baseline recruitment rate")

recruitment <- as.data.frame(postSum_list$R_year)
recruitment <- recruitment %>%
  select(area, Median, uCI, lCI, SD, CV) %>%
  mutate(Species = "Ptarmigan",
         Variable = "R",
         Definition = "Recruitment rate")

popdens <- as.data.frame(postSum_list$totDens_raw)
popdens <- popdens %>%
  select(area, Median, uCI, lCI, SD, CV) %>% 
  mutate(Median = Median * 1000000,
         uCI = uCI * 1000000,
         lCI = lCI * 1000000,
         SD = SD * 1000000,
         Species = "Ptarmigan",
         Variable = "totDens",
         Definition = "Population Density (birds/km^2)")

occ_bl <- as.data.frame(postSum_list$alphaPtar.Occ)
occ_bl <- occ_bl %>%
  select(area, Median, uCI, lCI, SD, CV) %>%
  mutate(Species = "Gyrfalcon",
         Variable = "μ_Init",
         Definition = "Baseline probability of brood initiation")

occ <- as.data.frame(postSum_list$probOcc)
occ <- occ %>%
  select(area, Median, uCI, lCI, SD, CV) %>%
  mutate(Species = "Gyrfalcon",
         Variable = "probInit",
         Definition = "Probability of brood initiation")

prod_bl <- as.data.frame(postSum_list$alphaPtar.Prod)
prod_bl <- prod_bl %>%
  select(area, Median, uCI, lCI, SD, CV) %>%
  mutate(Species = "Gyrfalcon",
         Variable = "μ_Prod",
         Definition = "Baseline no. nestlings per breeding attempt")

prod <- as.data.frame(postSum_list$terrProd)
prod <- prod %>%
  select(area, Median, uCI, lCI, SD, CV) %>%
  mutate(Species = "Gyrfalcon",
         Variable = "terrProd",
         Definition = "No. nestlings per breeding attempt")

table_summary <- rbind(survival_bl, survival, recruitment_bl, recruitment,
                       popdens, occ_bl, occ, prod_bl, prod)

table_summary <- table_summary %>%
  mutate(CI = paste0("[", round(lCI, 2), ", ", round(uCI, 2), "]"),
         Median = round(Median, 3),
         Area = recode(area, 
                       `1` = "Hardangervidda",
                       `2` = "Dovrefjell",
                       `3` = "Børgefjell"))

wide_table_summary <- table_summary %>%
  select(Species, Variable, Definition, Area, Median, CI) %>%
  pivot_wider(
    id_cols = c(Species, Variable, Definition),
    names_from  = Area,
    values_from = c(Median, CI),
    names_glue  = "{Area}_{.value}"   # area first, then stat
  )

# reorder to get the right headers
stats <- c("Median", "CI")
stat_order <- c("Median", "CI")
new_order <- c("Species", "Variable", "Definition",
               as.vector(t(outer(areas, stats, paste, sep = "_"))))
wide_table_summary <- wide_table_summary[, new_order]

col_keys <- names(wide_table_summary)
line1 <- ifelse(col_keys %in% c("Species","Variable", "Definition"), col_keys,
                sub("_.*$", "", col_keys))   # Area
line2 <- ifelse(col_keys %in% c("Species","Variable", "Definition"), "",
                sub(".*_", "", col_keys))    # Stat

header_map <- data.frame(
  col_keys = col_keys,
  line1 = line1,
  line2 = line2,
  stringsAsFactors = FALSE
)

ft <- flextable(wide_table_summary) %>%
  set_header_df(mapping = header_map, key = "col_keys") %>%
  merge_h(part = "header") %>%    # merges area names across stats
  merge_v(part = "header") %>%
  merge_v(j = "Species") %>%      # merge species names vertically
  theme_vanilla() %>%
  autofit() %>%
  align(i = 1, j = 3:ncol(wide_table_summary), align = "center", part = "header")  # top header

