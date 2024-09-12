#' Extracts and visualizes two measures of generation time for each area. 
#' 
#' The first measure of generation time is defined as the inverse of the 
#' asymptotic growth rate's elasticity to changes in fecundity (juvenile
#' production / recruitment) as proposed by Brooks & Lebreton (2001).
#' The second measure of generation time is the per-generation growth rate,
#' calculated using Rage::gen_time(). 
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param N_areas integer. Number of areas included in analyses. 
#' @param area_names character vector containing area/location names.
#' @param area_coord a tibble containing coordinate data for each area in analysis.
#' @param mapNM sf / dataframe object containing map of Norwegian municipalities. 
#' @param save logical. If TRUE (default) saves posterior summaries of generation
#' time estimates as an .rds file in the root directory. 
#'
#' @return a data frame containing posterior summaries of estimates of 
#' generation time. 
#' @export
#'
#' @examples
#' 
extract_GenTime <- function(mcmc.out, 
                            N_areas, area_names,
                            area_coord, mapNM,
                            save = TRUE){
  
  
  ## Convert posterior samples to matrix format
  out.mat <- as.matrix(mcmc.out)
  
  ## Set up dataframe for storing results
  GT_data_all <- data.frame()
  
  for(x in 1:N_areas){
    
    ## Calculate posterior for area-specific generation time using two different approaches
    GT_elasF <- GT_R0 <- rep(NA, nrow(out.mat))
    
    for(i in 1:nrow(out.mat)){
      
      # Extract average survival & recruitment
      S <- out.mat[i, paste0("Mu.S[", x, "]")]
      R <- out.mat[i, paste0("Mu.R[", x, "]")]
      
      # 2 age-class matrix
      mat_age2 <- matrix(NA, nrow = 2, ncol = 2)
      mat_age2[1, 1:2] <- S*R
      mat_age2[2, 1:2] <- S
      
      # Calculate elasticity
      mat_elas <- popbio::elasticity(mat_age2)
      
      # Calculate generation time (inverse of recruitment elasticity)
      GT_elasF[i] <- 1 / sum(mat_elas[1,])
      
      # Split matrix into growth-survival and reproduction sub-matrices
      matU <- matR <- mat_age2
      matU[1, ] <- 0
      matR[2, ] <- 0
      
      # Calculate generation time (per-generation growth rate)
      GT_R0[i] <- Rage::gen_time(matU = matU, matR = matR)
      
    }
    
    ## Summarize posterior
    GT_data_area <- data.frame(Area = area_names[x],
                               GT_measure = c("Inverse of fecundity elasticity", "Per generation growth rate"),
                               Median = c(median(GT_elasF), median(GT_R0)),
                               lCI = c(quantile(GT_elasF, probs = 0.025), quantile(GT_R0, probs = 0.025)),
                               uCI = c(quantile(GT_elasF, probs = 0.975), quantile(GT_R0, probs = 0.975)),
                               Mean = c(mean(GT_elasF), mean(GT_R0)),
                               SD = c(sd(GT_elasF), sd(GT_R0)),
                               CV = c(sd(GT_elasF / mean(GT_elasF)), sd(GT_elasF / mean(GT_elasF))))
    
    GT_data_all <- rbind(GT_data_all, GT_data_area)
    
  }
  
  ## Make plotting directory
  ifelse(!dir.exists("Plots/GenerationTime"), dir.create("Plots/GenerationTime"), FALSE)
  
  ## Coordinate plot
  area_coord <- area_coord %>%
    dplyr::rename(Area = spatialUnit,
                  Longitude = longitudeAvg,
                  Latitude = latitudeAvg)
  
  GT_data_all <- GT_data_all %>%
    dplyr::left_join(area_coord, by = "Area")
  
  p_latGT <- ggplot(GT_data_all) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0, fatten = 4, alpha = 0.5) +
    geom_point(aes(x = Latitude, y = Median, fill = Longitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    ylab("Estimate") + 
    ggtitle("Generation time") + 
    facet_wrap(~GT_measure, scales = "free_y") + 
    theme_classic()
    
  pdf("Plots/GenerationTime/GenerationTime_Latitude.pdf", width = 8, height = 4)
  print(p_latGT)
  dev.off()
  

  ## Map plots
  mapNM.GT_elasF <- mapNM %>%
    dplyr::left_join(., subset(GT_data_all, GT_measure == "Inverse of fecundity elasticity"), by = "Area") 
  
  mapNM.GT_R0 <- mapNM %>%
    dplyr::left_join(., subset(GT_data_all, GT_measure == "Per generation growth rate"), by = "Area") 
  
  pdf("Plots/GenerationTime/GenerationTime_elasF_Map.pdf", width = 5, height = 6)
  print(
    tmap::tm_shape(mapNM.GT_elasF) + tmap::tm_polygons("Median", palette = "plasma", style = "cont", colorNA = "grey80")
  )
  print(
    tmap::tm_shape(mapNM.GT_elasF) + tmap::tm_polygons("SD", palette = rev("BuGn"), style = "cont", colorNA = "grey80")
  )
  print(
    tmap::tm_shape(mapNM.GT_elasF) + tmap::tm_polygons("CV", palette = rev("BuGn"), style = "cont", colorNA = "grey80")
  )
  dev.off()
  
  pdf("Plots/GenerationTime/GenerationTime_R0_Map.pdf", width = 5, height = 6)
  print(
    tmap::tm_shape(mapNM.GT_R0) + tmap::tm_polygons("Median", palette = "plasma", style = "cont", colorNA = "grey80")
  )
  print(
    tmap::tm_shape(mapNM.GT_R0) + tmap::tm_polygons("SD", palette = rev("BuGn"), style = "cont", colorNA = "grey80")
  )
  print(
    tmap::tm_shape(mapNM.GT_R0) + tmap::tm_polygons("CV", palette = rev("BuGn"), style = "cont", colorNA = "grey80")
  )
  dev.off()
  
  ## Save (optional) and return
  if(save){
    saveRDS(GT_data_all, file = "PosteriorSummaries_GenTime_byArea.rds")
  }
  
  return(GT_data_all)
}
