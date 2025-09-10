#' Simulate values for initializing integrated model
#'
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#' @param R_perF logical. If TRUE, treats recruitment rate as juvenile per adult female.
#' If FALSE, treats recruitment rate as juvenile per adult (sum of both sexes).
#' @param survVarT logical. If TRUE, survival is simulated including annual variation.
#' @param fitRodentCov logical. If TRUE, initial values are generated for rodent 
#' covariate (effect) and covariate is included in data simulation.
#' 
#' @return A list containing one complete set of initial values for the model.
#' @export
#'
#' @examples

simulateInits_Integ <- function(nim.data, nim.constants, R_perF, survVarT, fitRodentCov, initVals.seed){
  
  set.seed(initVals.seed)
  
  # Limits and constants #
  #----------------------#
  
  N_areas <- nim.constants$N_areas
  N_ageC <- nim.constants$N_ageC
  N_years <- nim.constants$N_years
  
  if(N_areas == 1){
    N_sites <- nim.constants$N_sites[1]
  }else{
    N_sites <- nim.constants$N_sites
  }
  
  L <- nim.data$L
  W <- nim.constants$W
  pi <- 3.141593
  A <- nim.data$A
  
  
  # Missing covariate values #
  #--------------------------#
  
  if(!is.null(nim.data$RodentOcc)){
    RodentOcc <- nim.data$RodentOcc
  }else{
    RodentOcc <- matrix(0, nrow = N_areas, ncol = N_years)
  }
  
  if(NA %in% RodentOcc){
    RodentOcc[which(is.na(RodentOcc))] <- rnorm(length(which(is.na(RodentOcc))), 0, 1)
  }
  
  Inits_RodentOcc <- RodentOcc
  Inits_RodentOcc[which(!is.na(nim.data$RodentOcc))] <- NA
  
  
  # Vital rates #
  #-------------#
  
  ## Area-specific survival parameters
  h.Mu.S <- runif(1, 0.40, 0.45)
  h.sigma.S <- runif(1, 0.05, 0.2)
  Mu.S1 <- runif(1, 0.5, 0.7)
  
  #mu.S <- EnvStats::rnormTrunc(N_areas, qlogis(h.Mu.S), sd = h.sigma.S, max = qlogis(Mu.S1))
  mu.S <- EnvStats::rnormTrunc(N_areas, qlogis(h.Mu.S), sd = h.sigma.S)
  
  # sigmaT.S <- runif(1, 0.05, 0.2)
  sigmaR.S <- runif(1, 0.05, 0.2)
  
  Mu.S <- rep(NA, N_areas)
  S <-  matrix(NA, nrow = N_areas, ncol = N_years-1)
  #S1 <- S2 <- rep(NA, N_years)
  
  #eps.S1.prop <- runif(1, 0.3, 0.8)
  
  if(survVarT){
    # epsT.S <- rep(0, N_years-1)
    #epsT.S <- rnorm(N_years-1, 0, sigmaT.S)
    epsR.S <- matrix(0, nrow = N_areas, ncol = N_years-1)
    #epsR.S <- matrix(rnorm(N_areas*N_years, 0, sigmaR.S), nrow = N_areas, ncol = N_years)
    
  }else{
    # epsT.S <- rep(0, N_years-1)
    epsR.S <- matrix(0, nrow = N_areas, ncol = N_years-1)
  }
  
  #for(x in 1:N_areas){
  #  Mu.S[x] <- plogis(mu.S[x])
  #  S[x, 1:(N_years-1)] <- plogis(qlogis(Mu.S[x]) + epsT.S[1:(N_years-1)] + epsR.S[x, 1:(N_years-1)])
  #}
  
  for(x in 1:N_areas){
    Mu.S[x] <- plogis(mu.S[x])
    #S[x, 1:(N_years-1)] <- plogis(qlogis(Mu.S[x]) + epsT.S + epsR.S[x, ])
    S[x, 1:(N_years-1)] <- plogis(qlogis(Mu.S[x]) + epsR.S[x, ])
  }
  
  #S1[1:(N_years-1)] <- plogis(qlogis(Mu.S1) + eps.S1.prop*(epsT.S[1:(N_years-1)] + epsR.S[nim.constants$SurvAreaIdx, 1:(N_years-1)]))
  #S2[1:(N_years-1)] <- S[nim.constants$SurvAreaIdx, 1:(N_years-1)] / S1[1:(N_years-1)]
  
  ## Area-specific reproductive parameters
  h.Mu.R  <- runif(1, 1.5, 3)
  h.sigma.R <- runif(1, 0.05, 0.2)

  h.Mu.betaR.R <- runif(1, 0.01, 0.1)
  h.sigma.betaR.R <- runif(1, 0.05, 0.1)
  
  Mu.R <- rlnorm(N_areas, meanlog = log(h.Mu.R), sdlog =  h.sigma.R)
  
  # if(fitRodentCov){
  #   betaR.R <- rnorm(N_areas, mean = h.Mu.betaR.R, sd = h.sigma.betaR.R)
  # }else{
  #   betaR.R <- rep(0, N_areas) # FOr area specific estimates
  # }
  
  if(fitRodentCov){
    betaR.R <- rnorm(1, mean = h.Mu.betaR.R, sd = h.sigma.betaR.R)  # single slope
  }else{
    betaR.R <- 0
  }
  
  
  # sigmaT.R <- runif(1, 0.05, 0.2)
  sigmaR.R <- runif(1, 0.05, 0.2)
  
  R_year <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  # epsT.R <- rep(0, N_years)
  #epsT.R <- rnorm(N_year, 0, sigmaT.R)
  epsR.R <- matrix(0, nrow = N_areas, ncol = N_years)
  #epsR.R <- matrix(rnorm(N_areas*N_years, 0, sigmaR.R), nrow = N_areas, ncol = N_years)
  
  for(x in 1:N_areas){
    #R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R[x]*RodentOcc[x, 1:N_years] + epsT.R[1:N_years] + epsR.R[x, 1:N_years]) 
    R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R * RodentOcc[x, 1:N_years] + epsR.R[x, 1:N_years]) # Global slope and no random area effect
  }
  
  
  
  # Detection parameters #
  #----------------------#
  
  ## Area-specific detection parameters
  h.mu.dd <- runif(1, 3.5, 5.5)
  h.sigma.dd <- runif(1, 0.05, 0.2)
  
  #mu.dd <- rnorm(N_areas, h.mu.dd, sd = h.sigma.dd)
  mu.dd <- rep(h.mu.dd, N_areas)
  
  # sigmaT.dd <- runif(1, 0.05, 0.2)
  sigmaR.dd <- runif(1, 0.05, 0.2)
  
  sigma <- esw <- p <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  
  # epsT.dd <- rep(0, N_years)
  #epsT.dd <- rnorm(N_years, 0, sd = sigmaT.dd)
  epsR.dd <- matrix(0, nrow = N_areas, ncol = N_years)
  #epsR.dd <- matrix(rnorm(N_areas*N_years, 0, sigmaR.dd), nrow = N_areas, ncol = N_years)
  
  for(x in 1:N_areas){
    # sigma[x, 1:N_years] <- exp(mu.dd[x] + epsT.dd[1:N_years] + epsR.dd[x, 1:N_years])
    sigma[x, 1:N_years] <- exp(mu.dd[x] + epsR.dd[x, 1:N_years])
  }
  
  for(x in 1:N_areas){
    esw[x, 1:N_years] <- sqrt(pi * sigma[x, 1:N_years]^2 / 2) 
    p[x, 1:N_years] <- min(esw[x, 1:N_years], W) / W
  }
  
  
  # Population model #
  #------------------#
  
  ## Initial densities / population sizes
  Mu.D1 <- rep(NA, N_areas)
  sigma.D <- runif(N_areas, 0.1, 2)
  
  N_exp <- Density <- array(0, dim = c(N_areas, N_ageC, max(N_sites), N_years))
  
  for(x in 1:N_areas){
    
    D_x_sum <- nim.data$N_a_line_year[x,2,,] / (L[x,,]*W*2)
    D_data <- D_x_sum[which(!is.na(D_x_sum) & D_x_sum > 0)]
    Mu.D1[x] <- runif(1, quantile(D_data, 0.25), quantile(D_data, 0.75))  
    
    for(j in 1:N_sites[x]){
      
      Density[x, 2, j, 1] <- Mu.D1[x]
      
      if(R_perF){
        Density[x, 1, j, 1] <- (Density[x, 2, j, 1]/2)*R_year[x, 1] # Juveniles 
      }else{
        Density[x, 1, j, 1] <- Density[x, 2, j, 1]*R_year[x, 1] # Juveniles
      }
      
      
      lambda1 <- Density[x, 1, j, 1]*L[x, j, 1]*W*2
      lambda2 <- Density[x, 2, j, 1]*L[x, j, 1]*W*2
      
      if (!is.na(lambda1) && lambda1 > 0) {
        N_exp[x, 1, j, 1] <- extraDistr::rtpois(1, lambda = lambda1, a = 1)
      } else {
        N_exp[x, 1, j, 1] <- 1  # Fallback value if lambda is NA or negative
      }
      
      if (!is.na(lambda2) && lambda2 > 0) {
        N_exp[x, 2, j, 1] <- extraDistr::rtpois(1, lambda = lambda2, a = 1)
      } else {
        N_exp[x, 2, j, 1] <- 1 # Fallback value if lambda is NA or negative
      }
      # 
      # N_exp[x, 1, j, 1] <- extraDistr::rtpois(1, lambda = Density[x, 1, j, 1]*L[x, j, 1]*W*2, a = 1)
      # N_exp[x, 2, j, 1] <- extraDistr::rtpois(1, lambda = Density[x, 2, j, 1]*L[x, j, 1]*W*2, a = 1)
    }
  }
  
  ## Population projection over time
  for(x in 1:N_areas){
    for(j in 1:N_sites[x]){
      for(t in 2:N_years){
        
        Density[x, 2, j, t] <- sum(Density[x, 1:N_ageC, j, t-1])*S[x, t-1] # Adults
        
        if(R_perF){
          Density[x, 1, j, t] <- (Density[x, 2, j, t]/2)*R_year[x, t] # Juveniles 
        }else{
          Density[x, 1, j, t] <- Density[x, 2, j, t]*R_year[x, t] # Juveniles
        }
        
        N_exp[x, 1:N_ageC, j, t] <- Density[x, 1:N_ageC, j, t]*L[x, j, t]*W*2
      }
    }
  }
  
  ## Area-specific population size and density
  N_tot_exp <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  for(x in 1:N_areas){
    for (t in 1:N_years){
      N_tot_exp[x, t] <- sum(N_exp[x, 1, 1:N_sites[x], t] + N_exp[x, 2, 1:N_sites[x], t])    ## Summing up expected number of birds in covered area; 
    }
  }
  
  ## Area-, year-, and age-class specific density (for monitoring)
  meanDens <- array(NA, dim = c(N_areas, N_ageC, N_years))
  
  for(x in 1:N_areas){
    for(a in 1:N_ageC){
      for(t in 1:N_years){
        meanDens[x, a, t] <- mean(Density[x, a, 1:N_sites[x], t])
      }
    } 
  }
  
  # Gyrfalcon model #
  #-----------------#
  
  terrOcc <- nim.data$terrOcc
  terrMonitoredOcc <- nim.data$terrMonitoredOcc
  chicksTot <- nim.data$chicksTot
  terrMonitoredProd <- nim.data$terrMonitoredProd
  
  # Ptarmigan covariate slopes
  betaPtar.Occ <- rnorm(1, 0, 1)
  betaPtar.Prod <- rnorm(1, 0, 1)
  
  # Intercepts
  alphaPtar.Occ <- runif(N_areas, 0, 1)
  alphaPtar.Prod <- rnorm(N_areas, 0, 1)
  
  # Productivity
  # Set this up if terrprod is used as a latent variable in the ptarmigan model
  # terrProd <- matrix(NA, nrow = N_areas, ncol = N_years)
  # 
  # for (x in 1:N_areas) {
  #   for (t in 2:N_years) {
  #     if (!is.na(chicksTot[x, t]) && terrMonitoredProd[x, t] > 0) {
  #       terrProd[x, t] <- chicksTot[x, t] / terrMonitoredProd[x, t]
  #     } else {
  #       terrProd[x, t] <- runif(1, 1, 3)  # fallback to a reasonable range
  #     }
  #   }
  # }
  
  # Occupancy
  # Set this up when probOcc is used as a latent variable elsewhere, using the model formula to define probOcc and calculate starting values
  # probOcc <- matrix(NA, nrow = N_areas, ncol = N_years)
  # 
  # for (x in 1:N_areas) {
  #   for (t in 2:N_years) {
  #     # Use observed occupancy proportion if available
  #     if (!is.na(terrOcc[x, t]) && terrMonitoredOcc[x, t] > 0) {
  #       probOcc[x, t] <- terrOcc[x, t] / terrMonitoredOcc[x, t]
  #     } else {
  #       probOcc[x, t] <- runif(1, 0.3, 0.7)  # fallback to a reasonable range
  #     }
  #   }
  # }
  
  
  # Random effects
  epsT.Occ <- rnorm(N_years, 0, 0.5)
  epsT.Prod <- rnorm(N_years, 0, 0.5)
  
  sigmaT.Occ <- runif(1, 0.1, 1)
  sigmaT.Prod <- runif(1, 0.1, 1)
  
  
  
  # Assembly #
  #----------#
  
  InitVals <- list(
    #b = runif(1, 1, 50), 
    
    Mu.D1 = Mu.D1, 
    sigma.D = sigma.D,
    eps.D1 = matrix(0, nrow = nim.constants$N_areas, ncol = max(N_sites)),
    
    Mu.R = Mu.R,
    h.Mu.betaR.R = h.Mu.betaR.R, h.sigma.betaR.R = h.sigma.betaR.R,
    h.Mu.R = h.Mu.R, h.sigma.R = h.sigma.R,
    # sigmaT.R = sigmaT.R, 
    sigmaR.R = sigmaR.R,
    # epsT.R = epsT.R, 
    epsR.R = epsR.R,
    epsA.R =  log(Mu.R) - log(h.Mu.R),
    R_year = R_year,
    
    mu.dd = mu.dd,
    h.mu.dd = h.mu.dd, h.sigma.dd = h.sigma.dd,
    # sigmaT.dd = sigmaT.dd, 
    sigmaR.dd = sigmaR.dd,
    # epsT.dd = epsT.dd, 
    epsR.dd = epsR.dd,
    epsA.dd = mu.dd - h.mu.dd,
    sigma = sigma, sigma2 = sigma^2,
    esw = esw,
    p = p,
    
    h.Mu.S = h.Mu.S,
    h.sigma.S = h.sigma.S,
    mu.S = mu.S, Mu.S = Mu.S, 
    # sigmaT.S = sigmaT.S, 
    sigmaR.S = sigmaR.S,
    # epsT.S = epsT.S, 
    epsR.S = epsR.S,
    epsA.S = mu.S - logit(h.Mu.S),
    #Mu.S1 = Mu.S1,
    #eps.S1.prop = eps.S1.prop,
    #S1 = S1, S2 = S2, 
    S = S,
    
    Density = Density,
    meanDens = meanDens,
    N_exp = N_exp,
    N_tot_exp = N_tot_exp,
    
    betaPtar.Occ = betaPtar.Occ,
    betaPtar.Prod = betaPtar.Prod,
    alphaPtar.Occ = alphaPtar.Occ,
    alphaPtar.Prod = alphaPtar.Prod,
    
    #terrProd = terrProd, # add later
    #probOcc = probOcc, # add later
    
    epsT.Occ = epsT.Occ,
    epsT.Prod = epsT.Prod,
    sigmaT.Occ = sigmaT.Occ,
    sigmaT.Prod = sigmaT.Prod
  )
  
  if(fitRodentCov){
    InitVals$h.Mu.betaR.R <- h.Mu.betaR.R
    InitVals$h.sigma.betaR.R <- h.sigma.betaR.R
    InitVals$betaR.R <- betaR.R
    InitVals$epsA.betaR.R <- betaR.R - h.Mu.betaR.R
    InitVals$RodentOcc <- Inits_RodentOcc
  }
  
  InitVals$betaGyr.S <- 0
  #* CRN: Initialized at 0 for now, but this may need changing before full integration. 
  
  return(InitVals)
}
