mySeed <- 0
set.seed(mySeed)

N_years <- 5
meanDens <- 28
sdDens <- 17

#------------------#
# DUMMY MODEL CODE #
#------------------#

dummy.code <- nimble::nimbleCode({
  
  ## Species A population model
  for(t in 2:N_years){
    density_A[1, t] <- density_A[2, t] * R_A[t]
    # Original: Density[x, 1, j, t] <- Density[x, 2, j, t]*R_year[x, t]

    density_A[2, t] <- sum(density_A[1:2, t-1]) * S_A[t-1]
    # Original: Density[x, 2, j, t] <- sum(Density[x, 1:N_ageC, j, t-1])*S[x, t-1] 
  }
  
  for(t in 1:N_years){
    totDens_A[t] <- density_A[1, t] + density_A[2, t]
    totDens_A_std[t] <- (totDens_A[t] - meanDens)/sdDens
    # Original: 
    # meanDens[x, a, t] <- sum(Density[x, a, 1:N_sites[x], t]) / N_sites[x]
    # totDens_raw[x, t] <- meanDens[x, 1, t] + meanDens[x, 2, t]
    # totDens_std[x, t] <- max(min(-10, (totDens_raw[x, t] - totDens_meanCov[x]) / totDens_sdCov[x]), 10) 
  }

  ## Species A vital rate models
  for(t in 1:N_years){
    log(R_A[t]) <- log(Mu.R_A)
    # Original:  R_year[x, t] <- exp(log(Mu.R[x]) + epsR.R[x, t])
    
    logit(S_A[t]) <- logit(Mu.S_A) + beta.vrB*VR_B[t]
    # Original:
    # GyrPressure[x, t] <- terrProd[x, t]
    # logit(S[x, t]) <- logit(Mu.S[x]) + betaGyr.S*GyrPressure[x, t]
  }
  
  ## Species B vital rate models
  VR_B[1] <- Mu.VR_B
  
  for(t in 2:N_years){
    log(VR_B[t]) <- log(Mu.VR_B) + beta.densA*totDens_A_std[t-1] 
    # Original: log(terrProd[x, t]) <- log(alphaPtar.Prod[x]) + betaPtar.Prod * totDens_std[x, t-1] + epsT.Prod[t]
  }
  
  ## Priors
  density_A[1, 1] <- density_A[2, 1] * R_A[1]
  # Original: Density[x, 1, j, 1] <- Density[x, 2, j, 1]*R_year[x, 1]
  
  density_A[2, 1] ~ dpois(5)
  # Original: Density[x, 2, j, 1] <- exp(log(Mu.D1[x]) + eps.D1[x, j])
  
  Mu.R_A ~ dpois(1.5)
  Mu.S_A ~ dbeta(7, 3)
  
  beta.vrB ~ dnorm(mean = -0.1, sd = 0.02)
  
  Mu.VR_B ~ dpois(2)
  beta.densA ~ dnorm(mean = 0.1, sd = 0.02)
  
})


#--------------------------------#
# DUMMY INITIAL VALUE SIMULATION #
#--------------------------------#

dummy.initSim <- function(N_years){
  
  # Set up vectors
  density_A <- matrix(NA, nrow = 2, ncol = N_years)
  R_A <- S_A <- VR_B <- totDens_A <- totDens_A_std <- rep(NA, N_years)
    
  # Set constant/first-year values
  Mu.S_A <- runif(1, 0.3, 0.7)
  Mu.R_A <- rpois(1, 1.5)
  Mu.VR_B <- rpois(1, 2)
  
  density_A[2, 1] <- round(runif(1, 3, 8))
  density_A[1, 1] <- density_A[2, 1] * Mu.R_A #(= R_A[1])
  
  VR_B[1] <- Mu.VR_B
  
  beta.vrB <- runif(1, -0.08, -0.02)
  beta.densA <- runif(1, 0.01, 0.05)
  
  # Calculate year-specific values
  for(t in 1:N_years){
    
    if(t > 1){
      VR_B[t] <- exp(log(Mu.VR_B) + beta.densA*totDens_A_std[t-1])
    }
    
    S_A[t] <- plogis(qlogis(Mu.S_A) + beta.vrB*VR_B[t])
    R_A[t] <- Mu.R_A
    
    if(t > 1){
      density_A[2, t] <- sum(density_A[1:2, t-1]) * S_A[t-1]
      density_A[1, t] <- density_A[2, t] * R_A[t]
    }
    
    totDens_A[t] <- sum(density_A[1:2, t])
    totDens_A_std[t] <- (totDens_A[t] - meanDens)/sdDens
  }
  
  initList <- list(
    density_A = density_A,
    totDens_A = totDens_A,
    
    R_A = R_A,
    S_A = S_A,
    Mu.R_A = Mu.R_A,
    Mu.S_A = Mu.S_A,
    
    beta.vrB = beta.vrB,
    beta.densA = beta.densA,
    
    VR_B = VR_B,
    Mu.VR_B = Mu.VR_B
  )
}

#---------------------#
# DUMMY MCMC SETTINGS #
#---------------------#

## MCMC parameters
nchains <- 1
niter <- 1000
nburnin <- 200
nthin <- 1

## Parameters to monitor
params <- c("density_A", "totDens_A", "totDens_A_std",
            "R_A", "Mu.R_A", "S_A", "Mu.S_A",
            "VR_B", "Mu.VR_B",
            "beta.vrB", "beta.densA")

## Sample initial values
dummy.inits <- dummy.initSim(N_years = N_years)

#----------------#
# DUMMY TEST RUN #
#----------------#

## Run model
dummy.out <- nimbleMCMC(code = dummy.code,
                        data = list(), 
                        constants = list(N_years = N_years,
                                         meanDens = meanDens,
                                         sdDens = sdDens),
                        inits = dummy.inits, 
                        monitors = params,
                        nchains = nchains, 
                        niter = niter, 
                        nburnin = nburnin, 
                        thin = nthin, 
                        samplesAsCodaMCMC = TRUE, 
                        setSeed = mySeed)

## Extract average and sd for total density
meanDens <- median(apply(dummy.out[,c(paste0("totDens_A[", 1:N_years, "]"))], 1, mean))
# --> ~ 28

sdDens <- median(apply(dummy.out[,c(paste0("totDens_A[", 1:N_years, "]"))], 1, sd))
# --> 17
