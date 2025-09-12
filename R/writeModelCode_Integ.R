#' Write integrated distance sampling model code
#'
#' @param survVarT logical. If TRUE, writes code for a model including random
#' year variation in survival probability. If FALSE, assumed constant survival
#' probability across time. 
#' @param telemetryData logical. If TRUE, uses information from telemetry data
#' from Lierne. If FALSE, only line transect data is used. 
#' @return an R call object specifying the model structure for integrated 
#' distance sampling model. 
#' @param fullLoopPP logical. If TRUE, two way interactions between ptarmigan and
#' gyrfalcon are included.If FALSE, only effect of ptarmigan density on gyrfalcon
#' is modelled explicitly while gyrfalcon effect on ptarmigan is included via 
#' external covariate.
#' @export
#'
#' @examples

writeModelCode_Integ <- function(survVarT, telemetryData, fullLoopPP){
  
  IDSM.code <- nimble::nimbleCode({
    
    # N_areas = number of areas
    # N_sites[x] = number of sites in area x
    # N_ageC = number of age classes
    # N_years = number of years
    # N_sumR_obs[x] = number of data points in juvenile:adult ratio counts
    
    # N_exp[x, a, j, t] = Number of age class a individuals in site j of area x in year t
    # Density[x, a, j, t] = Density of age class a individuals in site j of area x in year t
    # L[x, j, t] = length of transect line in site j of area x in year t
    # W = truncation distance for line transect surveys
    
    # Mu.D1[x] = average initial density in area x
    
    # S[x, t] = annual survival from year t to t+1 in area x
    # R_year[x, t] = recruitment rate in year t in area x
    # p[x, t] = average distance sampling detection rate in area x in year t
    # sigma[x, t] = average distance sampling detection decay rate in area x in year t
    
    # eps.D1[x, j] = random site effect on initial density area x (site j)
    
    
    ####################
    # POPULATION MODEL #
    ####################
    
    for(x in 1:N_areas){
      
      #-----------------------------------------#
      # Initial population size/density (t = 1) #
      #-----------------------------------------#
      
      for (j in 1:N_sites[x]){
        
        ## Adult densities
        Density[x, 2, j, 1] <- exp(log(Mu.D1[x]) + eps.D1[x, j])
        
        ## Juvenile densities
        if(R_perF){
          Density[x, 1, j, 1] <- (Density[x, 2, j, 1]/2)*R_year[x, 1] 
        }else{
          Density[x, 1, j, 1] <- Density[x, 2, j, 1]*R_year[x, 1]
        }
        
        ## Adult and juvenile numbers
        N_exp[x, 1:N_ageC, j, 1] <- Density[x, 1:N_ageC, j, 1]*L[x, j, 1]*W*2 
        
      } # j
      
      #-------------------------------#
      # Population dynamics for t > 1 #
      #-------------------------------#
      
      for(j in 1:N_sites[x]){
        for(t in 2:N_years){
          
          ## Adult densities
          Density[x, 2, j, t] <- sum(Density[x, 1:N_ageC, j, t-1])*S[x, t-1] 
          
          ## Juvenile densities
          if(R_perF){
            Density[x, 1, j, t] <- (Density[x, 2, j, t]/2)*R_year[x, t] 
          }else{
            Density[x, 1, j, t] <- Density[x, 2, j, t]*R_year[x, t]
          }
          
          ## Adult and juvenile numbers
          N_exp[x, 1:N_ageC, j, t] <- Density[x, 1:N_ageC, j, t]*L[x, j, t]*W*2
          
        } # t
      } # j
      
      
      #--------------------#
      # Derived parameters #
      #--------------------#
      
      ## Area- and year-specific total densities (numbers)
      for (t in 1:N_years){
        N_tot_exp[x, t] <- sum(N_exp[x, 1, 1:N_sites[x], t] + N_exp[x, 2, 1:N_sites[x], t])
      } # t
      
      ## Area-, year-, and age-class specific density (for monitoring)
      for (a in 1:N_ageC) {
        for (t in 1:N_years) {
          meanDens[x, a, t] <- sum(Density[x, a, 1:N_sites[x], t]) / N_sites[x]
        } # t
      } # a
    } # x
    
    ## Area and year specific total densities (latent variable)
    for (x in 1:N_areas){
      for(t in 1:N_years){
        totDens_raw[x, t] <- meanDens[x, 1, t] + meanDens[x, 2, t]
        totDens_std[x, t] <- max(min(10, (totDens_raw[x, t] - totDens_meanCov[x]) / totDens_sdCov[x]), -10) # Standardized
      } # t
    } # x
    
    
    #---------------------#
    # Gyrfalcon models    #
    #---------------------#
    
    for (x in 1:N_areas){
      
      ## Year 1 (no ptarmigan density estimate available)
      # Occupancy
      logit(probOcc[x, 1]) <- logit(alphaPtar.Occ[x]) + epsT.Occ[1]
      
      # Productivity
      log(terrProd[x, 1]) <- log(alphaPtar.Prod[x]) + epsT.Prod[1]
      
      
      ## Years 2+ (ptarmigan density estimate available)
      
      for (t in 2:N_years){
        # Occupancy
        logit(probOcc[x, t]) <- logit(alphaPtar.Occ[x]) + betaPtar.Occ * totDens_std[x, t-1] + epsT.Occ[t]

        # Productivity
        log(terrProd[x, t]) <- log(alphaPtar.Prod[x]) + betaPtar.Prod * totDens_std[x, t-1] + epsT.Prod[t]

      } # t

    } # x
    
    ## Gyrfalcon pressure covariate
    if(fullLoopPP){
      for(x in 1:N_areas){
          
          # For year 1:
          GyrPressure[x, 1] <- 2 * probOcc[x, 1] * terrMonitoredOcc[x, 1] + terrProd[x, 1] # Avoid using the time lag for t=1 
          
          # For years 2+:
          for(t in 2:N_years){
          GyrPressure[x, t] <- 0.5 * (2 * probOcc[x, t-1] * terrMonitoredOcc[x, t-1]) + # Number of gyrfalcons present in the first half of the ptarmigan 'year'
                               0.5 * (2 * probOcc[x, t] * terrMonitoredOcc[x, t]) + terrProd[x, t] # Number of gyrfalcons present in second half of the ptarmigan 'year'
          
        }
      }
    }
    
    ####################
    # DATA LIKELIHOODS #
    ####################
    
    for(x in 1:N_areas){
      
      ## Age-specific line transect counts
      # N_a_line_year[x, a, j, t] = number of age class a individuals detected in site j of area x in year t
      for(j in 1:N_sites[x]){
        for(t in 1:N_years){
          for(a in 1:N_ageC){
            
            N_a_line_year[x, a, j, t] ~ dpois(p[x, t]*N_exp[x, a, j, t])
            
          }
        }
      }
      
      
      
      ## Juvenile:adult ratios from line transect observations
      # N_sumR_obs[x] = number of observations in juvenile:adult count data in area x
      # sumR_obs[x, i] = i'th entry in juvenile count data for area x
      # sumAd_obs[x, i] = i'th entry in adult count data for area x
      for (i in 1:N_sumR_obs[x]){
        
        sumR_obs[x, i] ~ dpois(R_year[x, sumR_obs_year[x, i]]*sumAd_obs[x, i])
      }
      
      
      
      ## Line transect observation distances (likelihood using nimbleDistance::dHN)
      # N_obs[x] = number of observations in detection distance data in area x
      # y[x, i] = i'th entry in detection distance data for area x
      for (i in 1:N_obs[x]){ 
        
        y[x, i] ~ dHN(sigma = sigma[x, Year_obs[x, i]], Xmax = W, point = 0)
      }
      
      
      ## Gyrfalcon models

      # Gyrfalcon occupancy (per area)
      
      for (t in 1:N_years){
        # terrOcc[x, t] = number of territories occupied per area
        # probOcc[x, t] = probability of occupancy
        # terrMonitoredOcc[x, t] = number of monitored territories for occupancy
        terrOcc[x, t] ~ dbin(prob = probOcc[x, t], size = terrMonitoredOcc[x, t]) 
        
      } # t
      
    } # x
    
    # Gyrfalcon productivity (per territory)
    # chicksObs[i] = chicks produced per territory
    # terrProd[x, t] = expected nr of chicks per territory
    # chickObs_area[i] = i'th entry of area index for a territory
    # chickObs_year[i] = i'th entry of year index for a territory
    
    for (i in 1:N_terr){
      chicksObs[i] ~ dpois(terrProd[chicksObs_area[i], chicksObs_year[i]])
    }
    
    
    ################################
    # PARAMETER MODELS/CONSTRAINTS #
    ################################
    
    for(x in 1:N_areas){
      
      ## Distance sampling detection parameters
      
      for(t in 1:N_years){
        
        # Detection decay
        log(sigma[x, t]) <- mu.dd[x]  + epsR.dd[x, t]
        
        sigma2[x, t] <- sigma[x, t] * sigma[x, t]
        
        # Effective strip width
        esw[x, t] <- sqrt(pi * sigma2[x, t] / 2) 
        
        # Average detection rate 
        p[x, t] <- min(esw[x, t], W) / W
      }
      
      
      
      ## Annual recruitment rates
      
      if(fitRodentCov){
        # R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R[x]*RodentOcc[x, 1:N_years] + epsR.R[x, 1:N_years]) # Area specific slopes
        R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R * RodentOcc[x, 1:N_years] + epsR.R[x, 1:N_years])
      }else{
        R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + epsR.R[x, 1:N_years]) 
      }
      
      
      ## Annual survival probabilities
      
      # if(survVarT){
      #   #logit(S[x, 1:(N_years-1)]) <- logit(Mu.S[x] + epsR.S[x, 1:(N_years-1)]) # Old version
      #   logit(S[x, 1:(N_years-1)]) <- logit(Mu.S[x]) + epsR.S[x, 1:(N_years-1)] + betaGyr.S[x]*GyrPressure[x, 1:(N_years-1)]
      # }else{
      #   logit(S[x, 1:(N_years-1)]) <- logit(Mu.S[x]) + betaGyr.S[x]*GyrPressure[x, 1:(N_years-1)]
      # } # Area specific effect of gyrpressure
      
      for(t in 1:(N_years-1)){
        if(survVarT){
          #logit(S[x, t]) <- logit(Mu.S[x] + epsR.S[x, t]) # Old version
          logit(S[x, t]) <- logit(Mu.S[x]) + betaGyr.S*GyrPressure[x, t] + epsR.S[x, t] 
        }else{
          logit(S[x, t]) <- logit(Mu.S[x]) + betaGyr.S*GyrPressure[x, t]
        }
      }
      
      # Experimenting with including snowdepth in survival estimates  
      #   logit(S[x, 1:(N_years-1)]) <- logit(Mu.S[x]) + epsR.S[x, 1:(N_years-1)] + betaGyr.S[x]*GyrPressure[x, 1:(N_years-1)] + betaSD.S[x] * SDPreBrood[x, 1:(N_years-1)]
      # }else{
      #   logit(S[x, 1:(N_years-1)]) <- logit(Mu.S[x]) + betaGyr.S[x]*GyrPressure[x, 1:(N_years-1)] + betaSD.S[x] * SDPreBrood[x, 1:(N_years-1)]
      # }
      
    } # x
    
    
    
    ###########
    # PRIORS  #
    ###########
    
    #-----------------------#
    # Intercepts / averages #
    #-----------------------#
    
    for(x in 1:N_areas){
      
      ## Initial density
      Mu.D1[x] ~ dunif(0, 10) # Original prior
      #Mu.D1[x] ~ dunif(0, 5)
      
      ## Recruitment fixed effects
      Mu.R[x] ~ dunif(0, 10) # Original prior
      #Mu.R[x] ~ dunif(0, 5) # Test
      #logMu.R[x] ~ dnorm(0.5, 1)
      #Mu.R[x] <- exp(logMu.R[x])
      
      
      ## Survival fixed effects
      Mu.S[x] ~ dunif(0, 1) # Original prior
      #logit.Mu.S[x] ~ dnorm(0, 1) # Test
      #logit.Mu.S[x] ~ dnorm(0, 0.5)
      #Mu.S[x] <- ilogit(logit.Mu.S[x])
      
      ## Detection fixed effects
      mu.dd[x] ~ dunif(-10, 100)
      #mu.dd[x] ~ dnorm(0, 2)
    }
    
    
    #----------------#
    # Random effects #
    #----------------#
    
    ## Standard deviations
    
    # Recruitment
    sigmaR.R ~ dunif(0, 5)
    #sigmaR.R ~ T(dnorm(0, 1), 0, )
    
    # Survival 
    if(survVarT){
      sigmaR.S ~ dunif(0, 5) # Original prior
      #sigmaR.S ~ dunif(0, 1)
    }
    
    # Detection
    sigmaR.dd ~ dunif(0, 20)
    
    # Initial density
    for(x in 1:N_areas){
      sigma.D[x] ~ dunif(0, 20)
    }
    
    # Gyrfalcon model
    sigmaT.Occ ~ dunif(0, 15)
    sigmaT.Prod ~ dunif(0, 15)
    
    ## Random effect levels
    
    # Shared year variation
    for (t in 1:N_years){
      # epsT.Gyr[t] ~ dnorm(0, sd = sigmaT.Gyr) # Shared random effect for gyr occ and prod
      epsT.Occ[t] ~ dnorm(0, sd = sigmaT.Occ)
      epsT.Prod[t] ~ dnorm(0, sd = sigmaT.Prod)
    }
    
    # Residual variation
    for(x in 1:N_areas){
      for (t in 1:N_years){
        
        epsR.R[x, t] ~ dnorm(0, sd = sigmaR.R) # Recruitment
        epsR.dd[x, t] ~ dnorm(0, sd = sigmaR.dd) # Detection
        
      }
    }
    
    for(x in 1:N_areas){
      for (t in 1:(N_years-1)){
        
        if(survVarT){
          epsR.S[x, t] ~ dnorm(0, sd = sigmaR.S) # Survival
        }
      }
    }
    
    # Site/transect variation
    for(x in 1:N_areas){
      for(j in 1:N_sites[x]){
        eps.D1[x, j] ~ dnorm(0, sd = sigma.D[x])
      }
    }
    
    #-------------------#
    # Covariate effects #
    #-------------------#
    
    ## Rodent effect on ptarmigan reproduction
    if(fitRodentCov){
      
      # for(x in 1:N_areas){
      #   betaR.R[x] ~ dunif(-5, 10) # Area specific slopes
      # }
      betaR.R ~ dunif(-5, 5)
    }
    
    ## Gyrfalcon effect on ptarmigan survival
    # for(x in 1:N_areas){
    #   betaGyr.S[x] ~ dunif(-10, 10) # Area specific slopes
    # } 
    betaGyr.S ~ dunif(-5, 5)

    ## Ptarmigan density effect on ptarmigan occupancy and productivity
    betaPtar.Prod ~ dunif(-5, 5)
    betaPtar.Occ ~ dunif(-5, 5)
    
    #-----------------#
    # Gyrfalcon model #
    #-----------------#
    
    for(x in 1:N_areas){
      alphaPtar.Occ[x] ~ dunif(0, 1)
      alphaPtar.Prod[x] ~ dunif(0, 8)
      #betaPtar.R[x] ~ dunif(-5, 5) # Area specific slope
    }


    #------------------#
    # Other parameters #
    #------------------#
    
    pi <- 3.141593
    
    
    ###############################
    # COVARIATE IMPUTATION MODELS #
    ###############################
    
    if(fitRodentCov){
      for(x in 1:N_areas){
        for (t in 1:N_years){
          
          RodentOcc[x, t] ~ dnorm(mean = 0, sd = 1)
        }
      }
    }
    
    
  })
  
  return(IDSM.code)
}