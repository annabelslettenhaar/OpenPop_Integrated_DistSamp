#' Calculate sampling correlation between (area-specific) survival and recruitment
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param N_areas integer. Number of areas included in analyses. 
#' @param area_names character vector containing area/location names.
#' @param area_coord a tibble containing coordinate data for each area in analysis.
#' @param min_years integer vector. Indices of first year with available data for each area/location. 
#' @param max_years integer vector. Indices of last year with available data for each area/location. 
#'
#' @return a vector of pdf file names. The plots and data can be found in Plots/VitalRate_corr and root directory, respectively.
#' @export
#'
#' @examples

checkVRcorrs <- function(mcmc.out, 
                         N_areas, area_names, area_coord,
                         min_years, max_years){
  
  ## Convert posterior samples to matrix format
  out.mat <- as.matrix(mcmc.out)
  
  ## Set up data frame for posterior summaries
  VR_corrs <- data.frame()
  
  ## Set up vectors for storing survival and recruitment values
  S_all <- R_all <- R_all2 <- c()
  
  ## Calculate posterior correlation coefficients - by area
  for(x in 1:N_areas){
    
    # Determine area-specific year range
    area_yearIdxs_S <- (min_years[x]:(max_years[x]-1))
    area_yearIdxs_R <- ((min_years[x]+1):max_years[x])
    
    # Extract posterior samples for area-specifc survival and recruitment
    S_vec <- as.vector(out.mat[, paste0("S[", x, ", ", area_yearIdxs_S, "]")])
    R_vec <- as.vector(out.mat[, paste0("R_year[", x, ", ", area_yearIdxs_R, "]")])
    R_vec2 <- as.vector(out.mat[, paste0("R_year[", x, ", ", area_yearIdxs_S, "]")])
    
    # Calculate and store area-specific correlation coefficients
    corResults <- data.frame(area = area_names[x], 
                             S1R2_cor = cor.test(S_vec, R_vec)$estimate,
                             S1R1_cor = cor.test(S_vec, R_vec2)$estimate)

    VR_corrs <- rbind(VR_corrs, corResults)
    
    # Add area-specific S and R to joint vectors
    S_all <- c(S_all, S_vec)
    R_all <- c(R_all, R_vec)
    R_all2 <- c(R_all2, R_vec2)
    
  }
  
  ## Calculate and store overall correlation coefficients
  corResults <- data.frame(area = "Overall", 
                           S1R2_cor = cor.test(S_all, R_all)$estimate,
                           S1R1_cor = cor.test(S_all, R_all2)$estimate)
  
  VR_corrs <- rbind(VR_corrs, corResults)
  
  ## Save results as RDS
  saveRDS(VR_corrs, file = "VR_corrCoef.rds")
  write.csv(VR_corrs, "VR_corrCoef.csv", row.names = FALSE)
  
  ## Merge in coordinate information
  colnames(area_coord) <- c("area", "Longitude", "Latitude")
  VR_corrs <- merge(VR_corrs, area_coord, by = "area", all.x = TRUE)
  
  ## Plot results
  p_VRcorr1 <- ggplot(subset(VR_corrs, area != "Overall")) + 
    geom_point(aes(x = Latitude, y = S1R2_cor, color = Longitude), size = 2, alpha = 0.75) +
    geom_hline(aes(yintercept = subset(VR_corrs, area == "Overall")$S1R2_cor), linetype = "solid", color = "purple") + 
    geom_hline(aes(yintercept = 0), linetype = "dotted") + 
    ylim(-1, 1) + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    xlab("Latitude") +
    ylab("Sample correlation") + 
    annotate(geom = "text", x = 64, y = 0.8, label = paste0("Overall coefficient: ", round(subset(VR_corrs, area == "Overall")$S1R2_cor, digits = 3)), color = "purple") + 
    theme_classic()
  
  p_VRcorr2 <- ggplot(subset(VR_corrs, area != "Overall")) + 
    geom_point(aes(x = Latitude, y = S1R1_cor, color = Longitude), size = 2, alpha = 0.75) +
    geom_hline(aes(yintercept = subset(VR_corrs, area == "Overall")$S1R1_cor), linetype = "solid", color = "purple") + 
    geom_hline(aes(yintercept = 0), linetype = "dotted") + 
    ylim(-1, 1) + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    xlab("Latitude") +
    ylab("Sample correlation") + 
    annotate(geom = "text", x = 64, y = 0.8, label = paste0("Overall coefficient: ", round(subset(VR_corrs, area == "Overall")$S1R1_cor, digits = 3)), color = "purple") + 
    theme_classic()
  
  ifelse(!dir.exists("Plots/VitalRate_corr"), dir.create("Plots/VitalRate_corr"), FALSE)
  
  pdf("Plots/VitalRate_corr/SurvRepCorr_Latitude.pdf", width = 8, height = 4)
  print(
    ggpubr::ggarrange(p_VRcorr1 + ggtitle("Survival[t] vs. Recruitment[t+1]") + theme(axis.title.x = element_blank()),
                      p_VRcorr2 + ggtitle("Survival[t] vs. Recruitment[t]") + theme(axis.title.x = element_blank()), 
                      heights = c(1, 1), ncol = 2, 
                      legend = "right", common.legend = TRUE)
  )
  dev.off()
  
  ## Return filepaths
  filepaths <- c("VR_corrCoef.rds", "VR_corrCoef.csv", "Plots/VitalRate_corr/SurvRepCorr_Latitude.pdf")
  return(filepaths)
}
