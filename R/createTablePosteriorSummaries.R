
createTablePosteriorSummaries <- function(postSum_list) {
  library(dplyr)
  library(tidyr)
  library(flextable)
  
  # Helper to process each component
  process_component <- function(data, species, variable, definition, scale = 1) {
    as.data.frame(data) %>%
      select(area, Median, uCI, lCI, SD, CV) %>%
      mutate(Median = Median * scale,
             uCI = uCI * scale,
             lCI = lCI * scale,
             SD = SD * scale,
             Species = species,
             Variable = variable,
             Definition = definition)
  }
  
  # Process all components
  survival_bl <- process_component(postSum_list$Mu.S, "Ptarmigan", "μ_S", "Baseline survival probability")
  survival <- process_component(postSum_list$S, "Ptarmigan", "S", "Survival probability")
  recruitment_bl <- process_component(postSum_list$Mu.R, "Ptarmigan", "μ_R", "Baseline recruitment rate")
  recruitment <- process_component(postSum_list$R_year, "Ptarmigan", "R", "Recruitment rate")
  popdens <- process_component(postSum_list$totDens_raw, "Ptarmigan", "totDens", "Population Density (birds/km^2)", scale = 1e6)
  occ_bl <- process_component(postSum_list$alphaPtar.Occ, "Gyrfalcon", "μ_Init", "Baseline probability of brood initiation")
  occ <- process_component(postSum_list$probOcc, "Gyrfalcon", "probInit", "Probability of brood initiation")
  prod_bl <- process_component(postSum_list$alphaPtar.Prod, "Gyrfalcon", "μ_Prod", "Baseline no. nestlings per breeding attempt")
  prod <- process_component(postSum_list$terrProd, "Gyrfalcon", "terrProd", "No. nestlings per breeding attempt")
  
  # Combine all
  table_summary <- bind_rows(survival_bl, survival, recruitment_bl, recruitment,
                             popdens, occ_bl, occ, prod_bl, prod) %>%
    mutate(CI = paste0("[", round(lCI, 2), ", ", round(uCI, 2), "]"),
           Median = round(Median, 3),
           Area = recode(area, `1` = "Hardangervidda", `2` = "Dovrefjell", `3` = "Børgefjell"))
  
  # Pivot wider
  wide_table_summary <- table_summary %>%
    select(Species, Variable, Definition, Area, Median, CI) %>%
    pivot_wider(
      id_cols = c(Species, Variable, Definition),
      names_from  = Area,
      values_from = c(Median, CI),
      names_glue  = "{Area}_{.value}"
    )
  
  # Reorder columns
  areas <- c("Hardangervidda", "Dovrefjell", "Børgefjell")
  stats <- c("Median", "CI")
  new_order <- c("Species", "Variable", "Definition",
                 as.vector(t(outer(areas, stats, paste, sep = "_"))))
  wide_table_summary <- wide_table_summary[, new_order]
  
  # Header mapping
  col_keys <- names(wide_table_summary)
  line1 <- ifelse(col_keys %in% c("Species","Variable", "Definition"), col_keys,
                  sub("_.*$", "", col_keys))
  line2 <- ifelse(col_keys %in% c("Species","Variable", "Definition"), "",
                  sub(".*_", "", col_keys))
  
  header_map <- data.frame(
    col_keys = col_keys,
    line1 = line1,
    line2 = line2,
    stringsAsFactors = FALSE
  )
  
  # Create flextable
  ft <- flextable(wide_table_summary) %>%
    set_header_df(mapping = header_map, key = "col_keys") %>%
    merge_h(part = "header") %>%
    merge_v(part = "header") %>%
    merge_v(j = "Species") %>%
    theme_vanilla() %>%
    autofit() %>%
    align(i = 1, j = 3:ncol(wide_table_summary), align = "center", part = "header")
  
  return(ft)
} 

  
