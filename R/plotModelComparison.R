#' Compare posterior densities of major parameters for two or more models
#'
#' @param modelPaths vector of strings. Paths to .rds files for each models' posterior samples
#' @param modelChars vector of strings. Characteristics of each model (identifier)
#' @param N_areas integer. Number of areas to plot comparison for. 
#' @param area_names character vector. Names of areas to plot comparison for.
#' @param N_sites integer vector. Number of sites per area.
#' @param N_years integer. Number of years to plot comparison for.
#' @param minYear integer. First year in the analysis. 
#' @param maxYear integer. Last year in the analysis. 
#' @param max_years integer vector. Last year of data in each area.
#' @param survAreaIdx integer. Area index for area(s) with seasonal telemetry data. 
#' @param plotPath string. Directory into which to save plots
#' @param returnData logical. If TRUE, return collated data from all models (default = FALSE)
#' 
#' @return
#' @export
#'
#' @examples

plotModelComparison <- function(modelPaths, modelChars, 
                                N_areas, area_names, N_sites, N_years, 
                                minYear, maxYear, max_years, survAreaIdx, 
                                plotPath, returnData = FALSE){

  if(length(modelPaths) != length(modelChars)){
    stop("Unequal number of model paths and model characteristics provided.")
  }
  
  ## Set number of models to compare
  nMod <- length(modelPaths)
  
  ## Make plotting directory if necessary
  dir.create(plotPath, showWarnings = FALSE)
  
  #--------------------#
  # AVERAGE PARAMETERS #
  #--------------------#
  
  ## List relevant parameters
  # Averages
  hyperParams <- c("h.Mu.R", "h.sigma.R", "sigmaT.R", "sigmaR.R",
                   "h.Mu.S", "h.sigma.S", "sigmaT.S", "sigmaR.S",
                   "h.mu.dd", "h.sigma.dd", "sigmaT.dd", "sigmaR.dd",
                   "h.Mu.betaR.R", "h.sigma.betaR.R", 
                   "Mu.S1", "Mu.S2")
  
  areaParams <- c(paste0("Mu.S[", 1:N_areas, "]"), 
                  paste0("Mu.R[", 1:N_areas, "]"),
                  paste0("mu.dd[", 1:N_areas, "]"),
                  paste0("betaR.R[", 1:N_areas, "]"))

  
  ## Collate samples of relevant data posterior data in a list
  data.list <- list()
  
  for(n in 1:nMod){
    mcmc.out <- readRDS(modelPaths[n])
    
    # Calculate second seasonal survival average
    for(i in 1:length(mcmc.out)){
      Mu.S2 <- mcmc.out[[i]][, paste0("Mu.S[", survAreaIdx, "]")]/mcmc.out[[i]][, "Mu.S1"]
      mcmc.out[[i]] <- as.mcmc(cbind(mcmc.out[[i]], Mu.S2))
    }
    
    # Convert posterior samples into matrix with only relevant parameters
    mcmc.mat <- as.matrix(mcmc.out)[,c(hyperParams, areaParams)]

    # Reformat samples
    out.data <- reshape2::melt(mcmc.mat)
    colnames(out.data) <- c("Sample", "Parameter", "Value")
    
    out.data$Model <- modelChars[n]
    data.list[[n]] <- out.data
  }
  
  ## Combine model data into single data frame
  data.all <- dplyr::bind_rows(data.list)
  
  ## Plot overlapping posterior densities for different parameter groups
  
  # Hyper parameters
  mains <- hyperParams[which(!(hyperParams %in% c("Mu.S1", "Mu.S2")))]
  
  pdf(paste0(plotPath, '/ModelComp_HyperParams.pdf'), width = 8, height = 5)
  print(
    ggplot(subset(data.all, Parameter %in% mains)) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
  dev.off()
  
  # Seasonal survival in Lierne
  seasonS <- c("Mu.S1", "Mu.S2")
  
  pdf(paste0(plotPath, '/ModelComp_SeasonSurv.pdf'), width = 8, height = 5)
  print(
    ggplot(subset(data.all, Parameter %in% seasonS & Value <= 1)) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = "free_y") +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
  dev.off()
  
  # Area-specific recruitment
  pdf(paste0(plotPath, '/ModelComp_AreaRecruitment.pdf'), width = 8*1.5, height = 5*1.5)
  print(
    ggplot(subset(data.all, Parameter %in% paste0("Mu.R[", 1:N_areas, "]"))) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
  dev.off()
  
  # Area-specific survival
  pdf(paste0(plotPath, '/ModelComp_AreaSurvival.pdf'), width = 8*1.5, height = 5*1.5)
  print(
    ggplot(subset(data.all, Parameter %in% paste0("Mu.S[", 1:N_areas, "]"))) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
  dev.off()
  
  # Area-specific detection
  pdf(paste0(plotPath, '/ModelComp_AreaDetection.pdf'), width = 8*1.5, height = 5*1.5)
  print(
    ggplot(subset(data.all, Parameter %in% paste0("mu.dd[", 1:N_areas, "]"))) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
  dev.off()
  
  # Area-specific rodent effect
  pdf(paste0(plotPath, '/ModelComp_AreaRodentEff.pdf'), width = 8*1.5, height = 5*1.5)
  print(
    ggplot(subset(data.all, Parameter %in% paste0("betaR.R[", 1:N_areas, "]"))) + 
      geom_density(aes(x = Value, color = Model, fill = Model), alpha = 0.5) + 
      facet_wrap(~Parameter, scales = 'free') +
      scale_fill_viridis(discrete = T) + 
      scale_color_viridis(discrete = T) + 
      theme_bw() + theme(panel.grid = element_blank())
  )
  dev.off()
  
  
  #------------------------#
  # TIME-SERIES PARAMETERS #
  #------------------------#
  
  ## Summarize posteriors for relevant parameters
  rRep <- pSurv <- data.frame()
  pDetect <- data.frame()
  popDens <- data.frame()
  
  for(n in 1:nMod){
    mcmc.out <- readRDS(modelPaths[n])
    out.mat <- as.matrix(mcmc.out)
    
    for(i in 1:N_areas){
      
      # Prepare matrices for temporary storage of results
      rRep.sum <- pSurv.sum <- data.frame()
      pDetect.sum <- data.frame()
      popDens.sum <-  data.frame()
      
      # Determine area-specific year range
      #area_yearIdxs <- (min_years[i]:max_years[i])
      area_yearIdxs <- (1:(maxYear-minYear+1))
      area_years <- area_yearIdxs + (minYear - 1)
      
      for(t in 1:length(area_years)){
        
        # Summarize annual reproductive rates
        rRep_name <- paste0("R_year[",  i, ", ", area_yearIdxs[t], "]")
        rRep_add <- data.frame(Model = modelChars[n],
                               Area = area_names[i],
                               Year = area_years[t], 
                               Median = median(out.mat[, rRep_name]),
                               lCI = unname(quantile(out.mat[, rRep_name], probs = 0.025)),
                               uCI = unname(quantile(out.mat[, rRep_name], probs = 0.975)))
        rRep.sum <- rbind(rRep.sum, rRep_add)
        
        # Summarize annual survival rates
        pSurv_name <- paste0("S[",  i, ", ", area_yearIdxs[t], "]")
        
        if(t < max_years[i]){
          pSurv_add <- data.frame(Model = modelChars[n],
                                  Area = area_names[i],
                                  Year = area_years[t], 
                                  Median = median(out.mat[, pSurv_name]),
                                  lCI = unname(quantile(out.mat[, pSurv_name], probs = 0.025)),
                                  uCI = unname(quantile(out.mat[, pSurv_name], probs = 0.975)))
        }else{
          pSurv_add <- data.frame(Model = modelChars[n],
                                  Area = area_names[i],
                                  Year = area_years[t], 
                                  Median = NA,
                                  lCI = NA,
                                  uCI = NA)
        }
        pSurv.sum <- rbind(pSurv.sum, pSurv_add)
        
        # Summarize annual detection probabilities
        pDetect_name <- paste0("p[",  i, ", ", area_yearIdxs[t], "]")
        pDetect_add <- data.frame(Model = modelChars[n],
                                  Area = area_names[i],
                                  Year = area_years[t], 
                                  Median = median(out.mat[, pDetect_name]),
                                  lCI = unname(quantile(out.mat[, pDetect_name], probs = 0.025)),
                                  uCI = unname(quantile(out.mat[, pDetect_name], probs = 0.975)))
        pDetect.sum <- rbind(pDetect.sum, pDetect_add)
        
        # Summarize annual average population densities
        popDens_juv <- out.mat[, paste0("meanDens[",  i, ", 1, ", area_yearIdxs[t], "]")]
        popDens_ad <- out.mat[, paste0("meanDens[",  i, ", 2, ", area_yearIdxs[t], "]")]

        popDens_mean <- popDens_juv + popDens_ad
        
        popDens_add <- data.frame(Model = modelChars[n],
                                  Area = area_names[i],
                                  Year = area_years[t], 
                                  Median = median(popDens_mean),
                                  lCI = unname(quantile(popDens_mean, probs = 0.025)),
                                  uCI = unname(quantile(popDens_mean, probs = 0.975)))
        popDens.sum <- rbind(popDens.sum, popDens_add)
      }  
      
      rRep <- rbind(rRep, rRep.sum)
      pSurv <- rbind(pSurv, pSurv.sum)
      pDetect <- rbind(pDetect, pDetect.sum)
      popDens <- rbind(popDens, popDens.sum) 
    }
  }

  ## Plot time series model comparison
  
  # Reproductive rates
  pdf(paste0(plotPath, "/ModelComp_tR.pdf"), width = 8, height = 5)
  for(i in 1:N_areas){
    
    p_rRep <- ggplot(subset(rRep, Area == area_names[i]), aes(x = Year, color = Model, fill = Model)) + 
      geom_line(aes(y = Median)) + 
      geom_ribbon(aes(ymin = lCI, ymax = uCI), alpha = 0.35, color = NA) +
      scale_x_continuous(breaks = c(minYear:maxYear), limits = c(minYear, maxYear)) + 
      #ylim(min(rRep$lCI), max(rRep$uCI)) + 
      ylab("Reproductive rate") +
      scale_color_viridis_d() + 
      scale_fill_viridis_d() + 
      ggtitle(area_names[i]) + 
      theme_bw() + 
      theme(panel.grid.minor = element_blank(), 
            axis.text.x = element_text(angle = 45, vjust = 0.75))
    
    print(p_rRep)
    
  }
  dev.off()
    
  # Survival probabilities
  pdf(paste0(plotPath, "/ModelComp_tS.pdf"), width = 8, height = 5)
  for(i in 1:N_areas){
    
    p_pSurv <- ggplot(subset(pSurv, Area == area_names[i]), aes(x = Year, color = Model, fill = Model)) + 
      geom_line(aes(y = Median)) + 
      geom_ribbon(aes(ymin = lCI, ymax = uCI), alpha = 0.35, color = NA) +
      scale_x_continuous(breaks = c(minYear:maxYear), limits = c(minYear, maxYear)) + 
      #ylim(min(rRep$lCI), max(rRep$uCI)) + 
      ylab("Survival probability") +
      scale_color_viridis_d() + 
      scale_fill_viridis_d() + 
      ggtitle(area_names[i]) + 
      theme_bw() + 
      theme(panel.grid.minor = element_blank(), 
            axis.text.x = element_text(angle = 45, vjust = 0.75))
    
    print(p_pSurv)
    
  }
  dev.off()
  
  # Detection probabilities
  pdf(paste0(plotPath, "/ModelComp_tDet.pdf"), width = 8, height = 5)
  for(i in 1:N_areas){
    
    p_pDetect <- ggplot(subset(pDetect, Area == area_names[i]), aes(x = Year, color = Model, fill = Model)) + 
      geom_line(aes(y = Median)) + 
      geom_ribbon(aes(ymin = lCI, ymax = uCI), alpha = 0.35, color = NA) +
      scale_x_continuous(breaks = c(minYear:maxYear), limits = c(minYear, maxYear)) + 
      #ylim(min(rRep$lCI), max(rRep$uCI)) + 
      ylab("Detection probability") +
      scale_color_viridis_d() + 
      scale_fill_viridis_d() + 
      ggtitle(area_names[i]) + 
      theme_bw() + 
      theme(panel.grid.minor = element_blank(), 
            axis.text.x = element_text(angle = 45, vjust = 0.75))
    
    print(p_pDetect)
    
  }
  dev.off()
  
  # Average population densities
  pdf(paste0(plotPath, "/ModelComp_tDens.pdf"), width = 8, height = 5)
  for(i in 1:N_areas){
    
    p_popDens <- ggplot(subset(popDens, Area == area_names[i]), aes(x = Year, color = Model, fill = Model)) + 
      geom_line(aes(y = Median)) + 
      geom_ribbon(aes(ymin = lCI, ymax = uCI), alpha = 0.35, color = NA) +
      scale_x_continuous(breaks = c(minYear:maxYear), limits = c(minYear, maxYear)) + 
      #ylim(min(rRep$lCI), max(rRep$uCI)) + 
      ylab(bquote("Average population density " (birds/km^2))) +
      scale_color_viridis_d() + 
      scale_fill_viridis_d() + 
      ggtitle(area_names[i]) + 
      theme_bw() + 
      theme(panel.grid.minor = element_blank(), 
            axis.text.x = element_text(angle = 45, vjust = 0.75))
    
    print(p_popDens)
    
  }
  dev.off()
  
  
  ## Return all data
  if(returnData){
    return(list(averages = data.all,
                rRep = rRep, 
                pSurv = pSurv, 
                pDetect = pDetect,
                popDens = popDens))
  }else{
    return()
  }

}
