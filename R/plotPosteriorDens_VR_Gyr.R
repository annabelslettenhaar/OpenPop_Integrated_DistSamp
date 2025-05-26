#' Plot posterior densities of vital rate parameters
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param N_areas integer. Number of areas included in analyses. 
#' @param area_names character vector containing area/location names.
#' @param N_years integer. Number of years included in analyses. 
#' @param minYear integer. The first year considered in analyses.
#' @param survAreaIdx integer. 
#' @param survVarT logical. If TRUE, plots time-dependent survival in addition 
#' to time-average survival.
#' @param fitRodentCov logical. If TRUE, interprets recruitment intercept as
#' baseline instead of average.
#' 
#' @return a vector of pdf plot names. The plots can be found in Plots/VR_PostDens.
#' @export
#'
#' @examples

plotPosteriorDens_VR_Gyr <- function(mcmc.out, N_areas, area_names, N_years, minYear, survVarT, fitRodentCov){

  ## List relevant parameters
  # Averages
  relParams <- c(paste0("Mu.S[", 1:N_areas, "]"), 
                 paste0("Mu.R[", 1:N_areas, "]"))
  
  # Time-dependents
  idxGrid <- expand.grid(c(1:N_areas), c(1:N_years))
  relParams_t <- paste0("R_year[", idxGrid[,1], ", ", idxGrid[,2], "]")
  
  if(survVarT){
    relParams_t <- c(relParams_t, paste0("S[", idxGrid[,1], ", ", idxGrid[,2], "]"))
  }
  
  ## Reduce posterior samples to relevant parameters only
  for(n in 1:length(mcmc.out)){
    mcmc.out[[n]] <- coda::as.mcmc(mcmc.out[[n]][,which(colnames(mcmc.out[[n]]) %in% c(relParams, relParams_t))])
  }

  ## Reformat samples
  out.data <- reshape2::melt(as.matrix(mcmc.out))
  colnames(out.data) <- c("Sample", "Parameter", "Value")
  
  ## Add area and year indices
  idx.data <- stringi::stri_extract_all_regex(out.data$Parameter, "[0-9]+", simplify = TRUE)
  colnames(idx.data) <- c("AreaIdx", "YearIdx")
  idx.data <- data.frame(idx.data)
  idx.data$AreaIdx <- as.numeric(idx.data$AreaIdx)
  idx.data$YearIdx <- as.numeric(idx.data$YearIdx)
  idx.data$Year <- idx.data$YearIdx + minYear - 1
  
  out.data <- cbind(out.data, idx.data)
  
  ## Add vital rate, season, and edited year labels
  out.data <- out.data %>%
    dplyr::mutate(VRcat = ifelse(grepl("R", Parameter, fixed = TRUE), "Recruitment", "Survival"),
                  Season = "Annual") %>%
    dplyr::mutate(plotYear = case_when(!is.na(YearIdx) ~ as.factor(Year),
                                       fitRodentCov & VRcat == "Recruitment" ~ "Baseline",
                                       TRUE ~ "Average"))

  ## Assemble partial plots for all areas
  p_Mu.S <- p_Mu.R <- p_R <- list()
  if(survVarT){p_S <- list()}
  
  for(i in 1:N_areas){
    # Average annual survival
    p_Mu.S[[i]] <- ggplot(subset(out.data, Parameter == paste0("Mu.S[", i, "]"))) + 
      geom_density(aes(x = Value), colour = viridis::mako(5)[3], fill = viridis::mako(5)[3], alpha = 0.5) + 
      xlab("") + ylab("Density") + 
      ggtitle("Average annual survival") + 
      theme_classic()
    
    # Average annual recruitment
    p_Mu.R[[i]] <- ggplot(subset(out.data, Parameter == paste0("Mu.R[", i, "]"))) + 
      geom_density(aes(x = Value), colour = "orange", fill = "orange", alpha = 0.5) + 
      xlab("") + ylab("Density") + 
      ggtitle(ifelse(fitRodentCov, "Baseline annual recruitment", "Average annual recruitment")) + 
      theme_classic()
    
    # Average + time-dependent annual survival
    if(survVarT){
      p_S[[i]] <- ggplot(subset(out.data, AreaIdx == i & VRcat == "Survival" & Season == "Annual")) + 
        geom_density(aes(x = Value, colour = plotYear, fill = plotYear), alpha = 0.5) + 
        scale_fill_manual(name = "Year", values = c(rev(colorspace::sequential_hcl(N_years-1, palette = "Peach")), viridis::mako(5)[3])) + 
        scale_color_manual(name = "Year", values = c(rev(colorspace::sequential_hcl(N_years-1, palette = "Peach")), viridis::mako(5)[3])) +
        xlab("") + ylab("Density") + 
        ggtitle("Survival probability") + 
        theme_classic()
    }
  
    # Average + time-dependent annual reproduction
    p_R[[i]] <- ggplot(subset(out.data, AreaIdx == i & VRcat == "Recruitment")) + 
      geom_density(aes(x = Value, colour = plotYear, fill = plotYear), alpha = 0.5) + 
      scale_fill_manual(name = "Year", values = c(rev(colorspace::sequential_hcl(N_years, palette = "Mint")), "orange")) + 
      scale_color_manual(name = "Year", values = c(rev(colorspace::sequential_hcl(N_years, palette = "Mint")), "orange")) +
      xlab("") + ylab("Density") + 
      ggtitle("Recruitment rate") + 
      theme_classic()
      
  }
  

  ## Save composite plots to pdf
  ifelse(!dir.exists("Plots/VR_PostDens"), dir.create("Plots/VR_PostDens"), FALSE) ## Check if folder exists, if not create folder
  
  # Simple average/baseline vital rates
  pdf("Plots/VR_PostDens/PostDens_Mu_S_R.pdf", height = 6, width = 7)
  for(i in 1:N_areas){
    
    title <- cowplot::ggdraw() + 
      cowplot::draw_label(area_names[[i]], fontface = 'bold', x = 0, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7))
    
      print(
        cowplot::plot_grid(
          title, cowplot::plot_grid(p_Mu.S[[i]], NULL, p_Mu.R[[i]], NULL, labels = c("A)", "B)", "C)", "")),
          ncol = 1,
          rel_heights = c(0.1, 1))
      )
    }

  
  dev.off()
  
  # Average survival + time-dependent recruitment
  pdf("Plots/VR_PostDens/PostDens_Mu_S_tR.pdf", height = 7, width = 8)
  for(i in 1:N_areas){
    
    title <- cowplot::ggdraw() + 
      cowplot::draw_label(area_names[[i]], fontface = 'bold', x = 0, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7))
    
      top_row <- cowplot::plot_grid(NULL, p_Mu.S[[i]], labels = c("A)", "B)"), label_size = 12)
    
    bottom_row <- p_R[[i]] + 
      theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank()) +
      guides(colour = guide_legend(nrow = 2, byrow = TRUE), fill = guide_legend(nrow = 2, byrow = TRUE))
    
    print(
      cowplot::plot_grid(
        title, cowplot::plot_grid(top_row, bottom_row, ncol = 1, labels = c("", "C)"), rel_heights = c(0.5, 0.7)),
        ncol = 1,
        rel_heights = c(0.1, 1))
    )
  }
  dev.off()
  
  # Time-dependent survival and recruitment
  if(survVarT){
    pdf("Plots/VR_PostDens/PostDens_tS_tR.pdf", height = 7, width = 8)
    for(i in 1:N_areas){
      
      title <- cowplot::ggdraw() + 
        cowplot::draw_label(area_names[[i]], fontface = 'bold', x = 0, hjust = 0) +
        theme(plot.margin = margin(0, 0, 0, 7))
      
      top_row <- p_S[[i]] + 
        guides(colour = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))
      bottom_row <- p_R[[i]] + 
        guides(colour = guide_legend(ncol = 2), fill = guide_legend(ncol = 2))
      
      print(
        cowplot::plot_grid(
          title, cowplot::plot_grid(top_row, bottom_row, ncol = 1, labels = c("A)", "B)")),
          ncol = 1,
          rel_heights = c(0.1, 1))
      )
    }
    dev.off()
  }

  
  ## Set and return plot path
  plot.paths <- c("Plots/VR_PostDens/PostDens_Mu_S_R.pdf", "Plots/VR_PostDens/PostDens_Mu_S_tR.pdf")
  if(survVarT){
    plot.paths <- c(plot.paths, "Plots/VR_PostDens/PostDens_tS_tR.pdf")
  }
  
  return(plot.paths)
}
