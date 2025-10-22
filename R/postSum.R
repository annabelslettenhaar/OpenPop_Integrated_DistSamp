
postSum <- function(samps, var_list, var_type, N_areas = NULL) {
  # samps: matrix or mcmc object (rows = posterior draws, cols = parameters)
  # var_list: vector of variable prefixes
  # var_type: vector of same length as var_list, one of "area_year", "area", "overall"
  # N_areas: number of areas (needed for "area_year" or "area" variables)
  
  if(!is.matrix(samps)) samps <- as.matrix(samps)
  
  results_list <- list()
  
  for(i in seq_along(var_list)) {
    var_prefix <- var_list[i]
    type <- var_type[i]
    
    if(type %in% c("area_year", "area") & is.null(N_areas)) {
      stop("N_areas must be provided for area/year or area variables")
    }
    
    # Initialize container for posterior draws for overall calculation
    overall_draws <- c()
    
    if(type == "area_year") {
      # per area × year: average across years per draw
      area_summaries <- lapply(1:N_areas, function(a) {
        cols_area <- grep(paste0("^", var_prefix, "\\[", a, ", "), colnames(samps), value = TRUE)
        mean_across <- rowMeans(samps[, cols_area, drop = FALSE])
        
        # store for overall summary
        overall_draws <<- cbind(overall_draws, mean_across)
        
        data.frame(
          area   = a,
          Median = median(mean_across),
          lCI    = quantile(mean_across, 0.025),
          uCI    = quantile(mean_across, 0.975),
          SD     = sd(mean_across),
          CV     = sd(mean_across) / median(mean_across)
        )
      })
      per_area_df <- do.call(rbind, area_summaries)
      
      # Overall summary across areas
      overall_mean_draws <- rowMeans(overall_draws)
      overall_summary <- data.frame(
        area = "overall",
        Median = median(overall_mean_draws),
        lCI    = quantile(overall_mean_draws, 0.025),
        uCI    = quantile(overall_mean_draws, 0.975),
        SD     = sd(overall_mean_draws),
        CV     = sd(overall_mean_draws) / median(overall_mean_draws)
      )
      
      results_list[[var_prefix]] <- rbind(per_area_df, overall_summary)
      
    } else if(type == "area") {
      # per area only: each column = area
      area_summaries <- lapply(1:N_areas, function(a) {
        cols_area <- grep(paste0("^", var_prefix, "\\[", a, "\\]"), colnames(samps), value = TRUE)
        draws <- samps[, cols_area, drop = FALSE]
        
        overall_draws <<- cbind(overall_draws, draws)
        
        data.frame(
          area   = a,
          Median = median(draws),
          lCI    = quantile(draws, 0.025),
          uCI    = quantile(draws, 0.975),
          SD     = sd(draws),
          CV     = sd(draws) / median(draws)
        )
      })
      per_area_df <- do.call(rbind, area_summaries)
      
      overall_mean_draws <- rowMeans(overall_draws)
      overall_summary <- data.frame(
        area = "overall",
        Median = median(overall_mean_draws),
        lCI    = quantile(overall_mean_draws, 0.025),
        uCI    = quantile(overall_mean_draws, 0.975),
        SD     = sd(overall_mean_draws),
        CV     = sd(overall_mean_draws) / median(overall_mean_draws)
      )
      
      results_list[[var_prefix]] <- rbind(per_area_df, overall_summary)
      
    } else if(type == "overall") {
      # single posterior per draw
      cols_var <- grep(paste0("^", var_prefix), colnames(samps), value = TRUE)
      draws <- samps[, cols_var, drop = FALSE]
      
      overall_summary <- data.frame(
        area = "overall",
        Median = median(draws),
        lCI    = quantile(draws, 0.025),
        uCI    = quantile(draws, 0.975),
        SD     = sd(draws),
        CV     = sd(draws) / median(draws)
      )
      
      results_list[[var_prefix]] <- overall_summary
    } else {
      stop(paste("Unknown type for variable", var_prefix))
    }
  }
  
  return(results_list)
}