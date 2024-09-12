#' Plot detection probability as a function of distance
#'
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param maxDist numeric. Maximum distance (in m) to plot detection for.
#' @param N_areas integer. Number of areas included in analyses. 
#' @param area_names character vector containing area/location names.
#'
#' @return a vector containing the pdf plot name. The plot can be found in Plots/DetectFunction.
#' @export
#'
#' @examples

plotDetectFunction <- function(mcmc.out,
                               maxDist,
                               N_areas, area_names){
  
  
  ## Convert posterior samples to matrix
  mcmc.mat <- as.matrix(mcmc.out)
  
  ## Make sequence of distances to predict for
  dist <- seq(0, maxDist, length.out = 100)

  ## Assemble dataframe for storing posterior summaries of predictions
  pred.data <- data.frame()
  
  for(i in 1:N_areas){
    
    ## Extract posterior samples of relevant parameters
    sigma_avg <- exp(mcmc.mat[, paste0("mu.dd[", i, "]")])

    ## Make, summarise, and store predictions for each covariate value
    for(x in 1:100){
      p.pred <- exp(-((dist[x]^2)/(2*sigma_avg^2)))
      pred.temp <- data.frame(Area = area_names[i], 
                              distance = dist[x], 
                              pred_Median = median(p.pred),
                              pred_lCI = unname(quantile(p.pred, probs = 0.025)),
                              pred_uCI = unname(quantile(p.pred, probs = 0.975)))
      pred.data <- rbind(pred.data, pred.temp)
    }
  }
  
  ## Plot predictions
  ifelse(!dir.exists("Plots/DetectFunction"), dir.create("Plots/DetectFunction"), FALSE) ## Check if folder exists, if not create folder
  
  pdf("Plots/DetectFunction/DetectionProb_distance.pdf", width = 6, height = 4) 
  for(i in 1:N_areas){
    print(
      ggplot(subset(pred.data, Area == area_names[i]), aes(x = distance, y = pred_Median)) +
        geom_line(color = "#856CEB") + 
        geom_ribbon(aes(ymin = pred_lCI, ymax = pred_uCI), alpha = 0.5, fill = "#856CEB") + 
        xlab("Distance from transect (m)") +
        ylab("Detection probability") + 
        ggtitle(area_names[i]) + 
        theme_classic()
    )
  }
  dev.off()
  
  ## Return plot path
  plot.paths <- "Plots/DetectFunction/DetectionProb_distance.pdf"
  return(plot.paths)
  
}