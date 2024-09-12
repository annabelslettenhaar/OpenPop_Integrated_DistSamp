#' Plot posterior summaries for vital rate and population parameters as a function of latitude
#'
#' @param PostSum.list a list containing 5 dataframes with summarised posteriors for several
#' parameters by area: average recruitment (rRep.sum), average survival 
#' (pSurv.sum), rodent effect slope (betaR.sum, = NA when fitRodentCov = FALSE), 
#' average population density (popDens.sum), and population growth rate 
#' (lambda.sum). 
#' @param area_coord a tibble containing coordinate data for each area in analysis.
#' @param minYear integer. First year in the analysis.
#' @param maxYear integer. Last year in the analysis.
#' @param fitRodentCov logical. If TRUE, makes plot for rodent effect slope in 
#' addition to other plots. 
#'
#' @return a vector of pdf plot names. The plots can be found in Plots/Latitude.
#' @export
#'
#' @examples

plotLatitude <- function(PostSum.list, 
                         area_coord,
                         minYear, maxYear,
                         fitRodentCov){
  
  ## Add coordinate information to summarized posterior data and plot
  area_coord <- area_coord %>%
    dplyr::rename(Area = spatialUnit,
                  Longitude = longitudeAvg,
                  Latitude = latitudeAvg)
  
  # Average detection parameters
  detect.sum <- PostSum.list$detect.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  p_detect <- ggplot(detect.sum) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0, fatten = 4, alpha = 0.5) +
    geom_point(aes(x = Latitude, y = Median, fill = Longitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  # Average reproductive rates
  rRep.sum <- PostSum.list$rRep.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  p_rRep <- ggplot(rRep.sum) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0, fatten = 4, alpha = 0.5) +
    geom_point(aes(x = Latitude, y = Median, fill = Longitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  # Average annual survival 
  pSurv.sum <- PostSum.list$pSurv.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  p_pSurv <- ggplot(pSurv.sum) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0, fatten = 4, alpha = 0.5) +
    geom_point(aes(x = Latitude, y = Median, fill = Longitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  ## Average annual survival vs. reproductive rates
  pSurv_rRep.sum <- pSurv.sum %>%
    dplyr::rename(surv.Median = Median,
                  surv.lCI = lCI,
                  surv.uCI = uCI) %>%
    dplyr::select(-Mean, -SD, -CV, -Latitude, -Longitude) %>%
    dplyr::left_join(rRep.sum, by = "Area")
  
  p_survRep1 <- ggplot(pSurv_rRep.sum) + 
    geom_pointrange(aes(x = surv.Median, y = Median, ymin = lCI, ymax = uCI, colour = Latitude), size = 0, fatten = 1, linewidth = 0.5, alpha = 0.375) +
    geom_pointrange(aes(x = surv.Median, y = Median, xmin = surv.lCI, xmax = surv.uCI, colour = Latitude), size = 0, fatten = 1, linewidth = 0.5, alpha = 0.375) +
    geom_point(aes(x = surv.Median, y = Median, fill = Latitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    xlab("Survival probability") +
    ylab("Recruitment rate") + 
    theme_classic()
  
  p_survRep2 <- ggplot(pSurv_rRep.sum) + 
    geom_point(aes(x = surv.Median, y = Median, colour = Latitude), size = 2, alpha = 0.75) +
    xlab("Survival probability") +
    ylab("Recruitment rate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  # Rodent effects 
  if(fitRodentCov){
    betaR.sum <- PostSum.list$betaR.sum %>%
      dplyr::left_join(area_coord, by = "Area")
    
    p_betaR <- ggplot(betaR.sum) + 
      geom_hline(aes(yintercept = 0), color = "grey70", linetype = "dotted") + 
      geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0, fatten = 4, alpha = 0.5) +
      geom_point(aes(x = Latitude, y = Median, fill = Longitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
      ylab("Estimate") + 
      paletteer::scale_color_paletteer_c("grDevices::Temps") + 
      paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
      theme_classic()
    
    p_betaR_long <- ggplot(betaR.sum) + 
      geom_hline(aes(yintercept = 0), color = "grey70", linetype = "dotted") + 
      geom_pointrange(aes(x = Longitude, y = Median, ymin = lCI, ymax = uCI, colour = Latitude), size = 0, fatten = 4, alpha = 0.5) +
      geom_point(aes(x = Longitude, y = Median, fill = Latitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
      ylab("Rodent effect on recruitment") + 
      paletteer::scale_color_paletteer_c("grDevices::Temps") + 
      paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
      theme_classic()
  }
  
  # Average population densities
  popDens.sum <- PostSum.list$popDens.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  p_popDens1 <- ggplot(subset(popDens.sum, SummaryPeriod == paste0(minYear, "-", maxYear))) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0, fatten = 4, alpha = 0.5) +
    geom_point(aes(x = Latitude, y = Median, fill = Longitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  p_popDens2 <- ggplot(subset(popDens.sum, SummaryPeriod != paste0(minYear, "-", maxYear))) + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0, fatten = 4, alpha = 0.5) +
    geom_point(aes(x = Latitude, y = Median, fill = Longitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  # Average population growth rates
  lambda.sum <- PostSum.list$lambda.sum %>%
    dplyr::left_join(area_coord, by = "Area")
  
  p_lambda1 <- ggplot(subset(lambda.sum, SummaryPeriod == paste0(minYear, "-", maxYear))) + 
    geom_hline(aes(yintercept = 1), color = "grey70", linetype = "dotted") + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0, fatten = 4, alpha = 0.5) +
    geom_point(aes(x = Latitude, y = Median, fill = Longitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  p_lambda2 <- ggplot(subset(lambda.sum, SummaryPeriod != paste0(minYear, "-", maxYear))) + 
    geom_hline(aes(yintercept = 1), color = "grey70", linetype = "dotted") + 
    geom_pointrange(aes(x = Latitude, y = Median, ymin = lCI, ymax = uCI, colour = Longitude), size = 0, fatten = 4, alpha = 0.5) +
    geom_point(aes(x = Latitude, y = Median, fill = Longitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    ylab("Estimate") + 
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    theme_classic()
  
  ## Average pop. densities vs. average pop. growth rates
  dens_lambda.sum <- popDens.sum %>%
    dplyr::rename(dens.Median = Median,
                  dens.lCI = lCI,
                  dens.uCI = uCI) %>%
    dplyr::select(-Mean, -SD, -CV, -Latitude, -Longitude) %>%
    dplyr::left_join(lambda.sum, by = c("Area", "SummaryPeriod"))
  
  p_lambdaDens1 <- ggplot(subset(dens_lambda.sum, SummaryPeriod == paste0(minYear, "-", maxYear))) + 
    geom_pointrange(aes(x = dens.Median, y = Median, ymin = lCI, ymax = uCI, colour = Latitude), size = 0, fatten = 1, linewidth = 0.5, alpha = 0.375) +
    geom_pointrange(aes(x = dens.Median, y = Median, xmin = dens.lCI, xmax = dens.uCI, colour = Latitude), size = 0, fatten = 1, linewidth = 0.5, alpha = 0.375) +
    geom_point(aes(x = dens.Median, y = Median, fill = Latitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    xlab("Population density") +
    ylab("Population growth rate") + 
    theme_classic()
  
  p_lambdaDens2 <- ggplot(subset(dens_lambda.sum, SummaryPeriod != paste0(minYear, "-", maxYear))) + 
    geom_pointrange(aes(x = dens.Median, y = Median, ymin = lCI, ymax = uCI, colour = Latitude), size = 0, fatten = 1, linewidth = 0.5, alpha = 0.375) +
    geom_pointrange(aes(x = dens.Median, y = Median, xmin = dens.lCI, xmax = dens.uCI, colour = Latitude), size = 0, fatten = 1, linewidth = 0.5, alpha = 0.375) +
    geom_point(aes(x = dens.Median, y = Median, fill = Latitude), shape = 21, color = "black", size = 2, alpha = 0.75) +
    paletteer::scale_color_paletteer_c("grDevices::Temps") + 
    paletteer::scale_fill_paletteer_c("grDevices::Temps") + 
    xlab("Population density") +
    ylab("Population growth rate") +  
    theme_classic()
  
  ## Make plotting directory if it does not exist yet
  ifelse(!dir.exists("Plots/Latitude"), dir.create("Plots/Latitude"), FALSE)
  
  ## Assemble plots (detection parameters)
  pdf("Plots/Latitude/DetectParams_Latitude.pdf", width = 8, height = 2)
  print(
    p_detect
  )
  dev.off()
  
  ## Assemble plots (vital rate parameters)
  if(fitRodentCov){
    pdf("Plots/Latitude/VitalRateParams_Latitude.pdf", width = 8, height = 6)
    print(
      ggpubr::ggarrange(p_pSurv + ggtitle("Average survival rate") + theme(axis.title.x = element_blank()),
                        p_rRep + ggtitle("Average recruitment rate") + theme(axis.title.x = element_blank()), 
                        p_betaR + ggtitle("Rodent effect on recruitment"),
                        heights = c(1, 1, 1.2), ncol = 1, 
                        legend = "right", common.legend = TRUE)
    )
    dev.off()
  }else{
    pdf("Plots/Latitude/VitalRateParams_Latitude.pdf", width = 8, height = 4)
    print(
      ggpubr::ggarrange(p_rRep + ggtitle("Average recruitment rate") + theme(axis.title.x = element_blank()), 
                        p_pSurv + ggtitle("Average survival rate"),  
                        heights = c(1, 1.2), ncol = 1, 
                        legend = "right", common.legend = TRUE)
    )
    dev.off()
  }

  ## Assemble plots (vital rate correlation)
  pdf("Plots/Latitude/SurvRep_Latitude.pdf", width = 5, height = 4)
  print(
    p_survRep1
  )
  print(
    p_survRep2
  )
  dev.off()
  
  ## Assemble plots (population-level parameters)
  pdf("Plots/Latitude/PopParams1_Latitude.pdf", width = 8, height = 4)
  print(
    ggpubr::ggarrange(p_popDens1 + ggtitle(paste0("Average population density (", minYear, "-", maxYear, ")")) + theme(axis.title.x = element_blank()),
                      p_lambda1 + ggtitle(paste0("Average population growth rate (", minYear, "-", maxYear, ")")),
                      heights = c(1, 1.1), ncol = 1, 
                      legend = "right", common.legend = TRUE)
  )
  dev.off()
  
  
  shared_sumPeriod <- subset(popDens.sum, SummaryPeriod != paste0(minYear, "-", maxYear))$SummaryPeriod[1]
  
  pdf("Plots/Latitude/PopParams2_Latitude.pdf", width = 8, height = 4)
  print(ggpubr::ggarrange(p_popDens2 + ggtitle(paste0("Average population density (", shared_sumPeriod, ")")) + theme(axis.title.x = element_blank()),
                          p_lambda2 + ggtitle(paste0("Average population growth rate (", shared_sumPeriod, ")")),
                          heights = c(1, 1.1), ncol = 1, 
                          legend = "right", common.legend = TRUE)
  )
  dev.off()
  
  ## Assemble plots (population-level parameter correlation)
  pdf("Plots/Latitude/LambdaDens1_Latitude.pdf", width = 5, height = 4)
  print(
    p_lambdaDens1
  )
  dev.off()
  
  pdf("Plots/Latitude/LambdaDens2_Latitude.pdf", width = 5, height = 4)
  print(
    p_lambdaDens2
  )
  dev.off()
  
  ## Assemble plots (manuscript array)
  if(fitRodentCov){
    pdf("Plots/Latitude/AllParams_Latitude.pdf", width = 8, height = 6)
    print(
      ggpubr::ggarrange(ggpubr::ggarrange(p_lambdaDens2 + theme(legend.position = "none") + ggtitle("A)"),
                                          p_survRep1 + theme(legend.position = "none") + ggtitle("B)"),
                                          nrow = 1), 
                        p_betaR_long + ggtitle("C)"),  
                        heights = c(1.5, 1), ncol = 1, 
                        legend = "right", common.legend = TRUE)
    )
    dev.off()
  }

  
  ## Write and return plot paths
  plot.paths <- c("Plots/Latitude/DetectParams_Latitude.pdf",
                  "Plots/Latitude/VitalRateParams_Latitude.pdf",
                  "Plots/Latitude/SurvRep_Latitude.pdf",
                  "Plots/Latitude/PopParams1_Latitude.pdf",
                  "Plots/Latitude/PopParams2_Latitude.pdf",
                  "Plots/Latitude/LambdaDens1_Latitude.pdf",
                  "Plots/Latitude/LambdaDens2_Latitude.pdf")
  
  if(fitRodentCov){
    plot.paths <- c(plot.paths, "Plots/Latitude/AllParams_Latitude.pdf")
  }
  
  return(plot.paths)
  
}