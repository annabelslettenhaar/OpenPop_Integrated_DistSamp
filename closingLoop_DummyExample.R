mySeed <- 0
set.seed(mySeed)

N_years <- 5


#------------------#
# DUMMY MODEL CODE #
#------------------#

dummy.code <- nimble::nimbleCode({
  
  ## Species A population model
  for(t in 2:N_years){
    density_A[t] <- density_A[t-1] * lambda_A[t-1]
  }
  
  ## Species A population growth rate
  for(t in 1:(N_years-1)){
    log(lambda_A[t]) <- log(Mu.lambda) + beta.vrB*VR_B[t]
  }
  
  ## Species B vital rate models
  VR_B[1] <- Mu.VR
  
  for(t in 2:N_years){
    log(VR_B[t]) <- log(Mu.VR) + beta.densA*density_A[t-1] 
  }
  
  ## Priors
  density_A[1] ~ dpois(3)

  Mu.lambda ~ dlnorm(meanlog = 0, sdlog = 1)
  beta.vrB ~ dnorm(mean = -0.1, sd = 0.02)
  
  Mu.VR ~ dpois(2)
  beta.densA ~ dnorm(mean = 0.1, sd = 0.02)
  
})


#--------------------------------#
# DUMMY INITIAL VALUE SIMULATION #
#--------------------------------#

dummy.initSim <- function(N_years){
  
  # Set up vectors
  density_A <- VR_B <- rep(NA, N_years)
  lambda_A <- rep(NA, N_years-1)
    
  # Set constant/first-year values
  Mu.lambda <- rlnorm(1, meanlog = 0, sdlog = 0.1)
  Mu.VR <- rpois(1, 2)
  
  density_A[1] <- round(runif(1, 2, 4))
  VR_B[1] <- Mu.VR
  
  beta.vrB <- runif(1, -0.08, -0.02)
  beta.densA <- runif(1, 0.01, 0.05)
  
  # Calculate year-specific values
  for(t in 2:N_years){
    
    VR_B[t] <- exp(log(Mu.VR) + beta.densA*density_A[t-1])
    lambda_A[t-1] <- exp(log(Mu.lambda) + beta.vrB*VR_B[t-1])
    
    density_A[t] <- density_A[t-1] * lambda_A[t-1]
  }
  
  initList <- list(
    density_A = density_A,
    
    lambda_A = lambda_A,
    Mu.lambda = Mu.lambda,
    
    beta.vrB = beta.vrB,
    beta.densA = beta.densA,
    
    VR_B = VR_B,
    Mu.VR = Mu.VR
  )
}

#---------------------#
# DUMMY MCMC SETTINGS #
#---------------------#

## MCMC parameters
nchains <- 1
niter <- 10
nburnin <- 0
nthin <- 1

## Parameters to monitor
params <- c("density_A", "lambda_A", "Mu.lambda", "VR_B", "Mu.VR")

## Sample initial values
dummy.inits <- dummy.initSim(N_years = N_years)

dummy.out <- nimbleMCMC(code = dummy.code,
                        data = list(), 
                        constants = list(N_years = N_years),
                        inits = dummy.inits, 
                        monitors = params,
                        nchains = nchains, 
                        niter = niter, 
                        nburnin = nburnin, 
                        thin = nthin, 
                        samplesAsCodaMCMC = TRUE, 
                        setSeed = mySeed)

