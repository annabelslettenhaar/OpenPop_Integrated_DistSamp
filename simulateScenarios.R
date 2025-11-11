
# Limits and constants #
#----------------------#
  
N_areas <- input_data$nim.constants$N_areas
N_ageC <- input_data$nim.constants$N_ageC
N_years <- input_data$nim.constants$N_years

if(N_areas == 1){
  N_sites <- input_data$nim.constants$N_sites[1]
}else{
  N_sites <- input_data$nim.constants$N_sites
}

L <- input_data$nim.data$L
W <- input_data$nim.constants$W
pi <- 3.141593
A <- input_data$nim.data$A

totDens_meanCov <- input_data$nim.constants$totDens_meanCov
totDens_sdCov <- input_data$nim.constants$totDens_sdCov


# Rodent data #
#-------------#

# Second order auto regressive model (AR2) to simulate rodent dynamics 

simulateRodentAR2 <- function(N_areas, N_years, alpha, phi1, phi2, sigma) {
  mat <- matrix(NA, nrow = N_areas, ncol = N_years)
  for (i in 1:N_areas) {
    mat[i, 1:2] <- rnorm(2, 0, 1)
    for (t in 3:N_years) {
      mat[i, t] <- alpha + phi1 * mat[i, t-1] + phi2 * mat[i, t-2] + rnorm(1, 0, sigma)
    }
  }
  mat
}

# Parameters with three versions of rodent fluctuation strength
params <- list(
  weak = list(alpha = 0, phi1 = -0.1, phi2 = -0.15, sigma = 0.5), # weak fluctuations
  moderate = list(alpha=0, phi1=0.6, phi2=-0.3, sigma=1.0), # moderate fluctuations
  strong = list(alpha=0, phi1=1.2, phi2=-0.5, sigma=1.5) # strong fluctuations
)

# Simulate and plot
par(mfrow=c(3,1))
for (scenario in names(params)) {
  p <- params[[scenario]]
  sim <- simulateRodentAR2(3, 30, p$alpha, p$phi1, p$phi2, p$sigma)
  matplot(t(sim), type='l', lty=1, main=paste("Scenario:", scenario),
          ylab="Rodent index", xlab="Year")
}

# Choose scenario and run simulation
chosen <- "strong"
p <- params[[chosen]]

# Run simulation
RodentOcc <- simulateRodentAR2(N_areas = 3, N_years = 30, 
                         alpha = p$alpha, phi1 = p$phi1, phi2 = p$phi2, sigma = p$sigma)


matplot(t(RodentOcc), type='l', lty=1, main=paste("Scenario:", chosen),
        ylab="Rodent index", xlab="Year")


# Gyrfalcon vital rates #
#-----------------------#

## Extract posterior medians
post_gyr <- extractPostMedians(modelOutput = model_output,
                               paramNames = c("alphaPtar.Occ", "alphaPtar.Prod", 
                                              "betaPtar.Occ", "betaPtar.Prod", 
                                              "sigmaT.Occ", "sigmaT.Prod"))

# Intercepts
alphaPtar.Occ <- post_gyr$alphaPtar.Occ 
alphaPtar.Prod <- post_gyr$alphaPtar.Prod 

# Ptarmigan covariate slopes
betaPtar.Occ <- post_gyr$betaPtar.Occ 
betaPtar.Prod <- post_gyr$betaPtar.Prod 

# Random effects
sigmaT.Occ <- post_gyr$sigmaT.Occ 
sigmaT.Prod <- post_gyr$sigmaT.Prod 

epsT.Occ <- rnorm(N_years, 0, sigmaT.Occ)
epsT.Prod <- rnorm(N_years, 0, sigmaT.Prod)

# Storing estimates
probOcc <- terrProd <- matrix(NA, nrow = N_areas, ncol = N_years)

# Initial values for year 1 
probOcc[, 1] <- plogis(alphaPtar.Occ)
terrProd[, 1] <- exp(alphaPtar.Prod)

# Replace this with the full model after everything is initialised
# for(x in 1:N_areas){
#   for(t in 1:N_years){
#     probOcc[x, t] <- plogis(qlogis(alphaPtar.Occ[x]) + epsT.Occ[t])
#     terrProd[x, t] <- exp(log(alphaPtar.Prod[x]) + epsT.Prod[t])
#   }
# }


# Ptarmigan vital rates #
#-----------------------#

## Extract posterior medians
post_ptar <- extractPostMedians(modelOutput = model_output,
                                paramNames = c("Mu.S", "Mu.R", 
                                               "sigmaR.S", "sigmaR.R", 
                                               "betaR.R", "betaTemp.R",
                                               "betaGyr.S"))
# Intercepts
Mu.S <- post_ptar$Mu.S
Mu.R <- post_ptar$Mu.R

# Covariate slopes
betaGyr.S <- post_ptar$betaGyr.S
betaTemp.R <- post_ptar$betaTemp.R

if(fitRodentCov){
  betaR.R <- post_ptar$betaR.R
}else{
  betaR.R <- 0
}

# Random effects
sigmaR.S <- post_ptar$sigmaR.S
sigmaR.R <- post_ptar$sigmaR.R

epsR.S <- matrix(rnorm(N_areas*(N_years-1), 0, sigmaR.S), nrow = N_areas)
epsR.R <- matrix(rnorm(N_areas*N_years, 0, sigmaR.R), nrow = N_areas)


# if(survVarT){
#   epsR.S <- matrix(0, nrow = N_areas, ncol = N_years-1)
#   #epsR.S <- matrix(rnorm(N_areas*N_years, 0, sigmaR.S), nrow = N_areas, ncol = N_years)
# }else{
#   epsR.S <- matrix(0, nrow = N_areas, ncol = N_years-1)
# }

# Storing and model survival
S <-  matrix(NA, nrow = N_areas, ncol = N_years-1)

for(x in 1:N_areas){
  S[x, 1:(N_years-1)] <- plogis(qlogis(Mu.S[x]) + epsR.S[x, ])
}

# Storing and model recruitment
R_year <- matrix(NA, nrow = N_areas, ncol = N_years)

for(x in 1:N_areas){
  R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R * RodentOcc[x, 1:N_years] + epsR.R[x, 1:N_years])
}


# Detection parameters #
#----------------------#

# Leave this whole block out and use model estimates instead

# Something like this?
ptarDens <- posterior_samples$ptardens


# Population model #
#------------------#

## Extract posterior median
post_pop <- extractPostMedians(modelOutput = model_output,
                              paramNames = c("Mu.D1"))


## Initial densities / population sizes

# Intercept
Mu.D1 <- post_pop$Mu.D1



# Putting models together #
#-------------------------#


## Gyrfalcon occupancy and productivity

for (t in 2:N_years) {
  for (x in 1:N_areas) {
    probOcc[x, t] <- plogis(alphaPtar.Occ[x] + betaPtar.Occ * totDens_std[x, t-1] + epsT.Occ[t])
    terrProd[x, t] <- exp(alphaPtar.Prod[x] + betaPtar.Prod * totDens_std[x, t-1] + epsT.Prod[t])
  }
}


## Gyrfalcon pressure covariate

for(x in 1:N_areas){
  
  # For year 1:
  GyrPressure_raw[x, 1] <- probOcc[x, 1] 
  
  # For years 2+:
  for(t in 2:N_years){
    GyrPressure_raw[x, t] <- (0.5 * probOcc[x, t-1]) + # Occupancy probability in the first half of the ptarmigan 'year'
      (0.5 * probOcc[x, t]) # Occupancy probability in second half of the ptarmigan 'year'
  }
  
  for(t in 1:N_years){
    GyrPressure_std[x, t] <- (GyrPressure_raw[x, t] - GyrPressure_meanCov[x]) / GyrPressure_sdCov[x] # Standardizing GyrPressure
  }
}


## Ptarmigan survival

for(x in 1:N_areas){
  S[x, 1:(N_years-1)] <- plogis(qlogis(Mu.S[x]) + betaGyr.S * GyrPressure_std[x, t] + 
                                  epsR.S[x, ])
}


## Ptarmigan recruitment

for(x in 1:N_areas){
  R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R * RodentOcc[x, 1:N_years] + 
                                betaTemp.R * SpringTemp[x, 1:N_years] + 
                                epsR.R[x, 1:N_years])
}

## Ptarmigan density

# Setup matrix
N_exp <- Density <- array(0, dim = c(N_areas, N_ageC, max(N_sites), N_years))

for(x in 1:N_areas){
  
  # D_x_sum <- nim.data$N_a_line_year[x,2,,] / (L[x,,]*W*2)
  # D_data <- D_x_sum[which(!is.na(D_x_sum) & D_x_sum > 0)]
  # Mu.D1[x] <- runif(1, quantile(D_data, 0.25), quantile(D_data, 0.75))  
  
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

## Density covariate
totDens_raw <- totDens_std <- matrix(NA, nrow = N_areas, ncol = N_years)

for (x in 1:N_areas){
  for(t in 1:N_years){
    totDens_raw[x, t] <- meanDens[x, 1, t] + meanDens[x, 2, t]
    totDens_std[x, t] <- max(min(-10, (totDens_raw[x, t] - totDens_meanCov[x]) / totDens_sdCov[x]), 10) # Standardized
  }
}


