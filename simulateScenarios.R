
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
RodentOcc <- simulateRodentAR2(N_areas = 3, N_years = 50, 
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

# # Storing estimates
# probOcc <- terrProd <- matrix(NA, nrow = N_areas, ncol = N_years)
# 
# # Initial values for year 1 
# probOcc[, 1] <- plogis(alphaPtar.Occ)
# terrProd[, 1] <- exp(alphaPtar.Prod)


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


# Population density #
#------------------#

## Extract posterior median
post_pop <- extractPostMedians(modelOutput = model_output,
                              paramNames = c("Mu.D1"))


## Initial densities / population sizes

# Intercept
Mu.D1 <- post_pop$Mu.D1


# Initialize matrices #
#---------------------#

# GyrPressure_raw <- GyrPressure_std <- matrix(NA, nrow = N_areas, ncol = N_years)
# S <- matrix(NA, nrow = N_areas, ncol = N_years-1)
# R_year <- matrix(NA, nrow = N_areas, ncol = N_years)
# Density <- array(0, dim = c(N_areas, N_ageC, max(N_sites), N_years))
# N_exp <- array(0, dim = c(N_areas, N_ageC, max(N_sites), N_years))
# meanDens <- array(NA, dim = c(N_areas, N_ageC, N_years))
# totDens_raw <- totDens_std <- matrix(NA, nrow = N_areas, ncol = N_years)


# Putting models together #
#-------------------------#

# Starting with the simplest of the simplest: only ptarmigan dynamics without the site level variation 
# No standardization and no random effects

alphaPtar.Occ <- -5
betaPtar.Occ <- 1
Mu.S <- 0.5
Mu.R <- 1.5
betaGyr.S <- -2

sim.years <- 50
N_years <- sim.years

# Initialize matrices
AdultDensity <- JuvenileDensity <- totalDensity <- matrix(NA, nrow = N_areas, ncol = N_years)
S <- matrix(NA, nrow = N_areas, ncol = N_years)
R_year <- matrix(NA, nrow = N_areas, ncol = N_years)
probOcc <- terrProd <- matrix(NA, nrow = N_areas, ncol = N_years)
GyrPressure <- matrix(NA, nrow = N_areas, ncol = N_years)

# Initial values for the first year
# AdultDensity[, 1] <- Mu.D1 * 1000000
AdultDensity[, 1] <- 10
# JuvenileDensity[, 1] <- if (R_perF) (AdultDensity[, 1]/2)*Mu.R else AdultDensity[, 1]*Mu.R
JuvenileDensity[, 1] <- 5
totalDensity[, 1] <- AdultDensity[, 1] + JuvenileDensity[, 1]
# probOcc[, 1] <- plogis(alphaPtar.Occ)
probOcc[, 1] <- 0.5
terrProd[, 1] <- exp(alphaPtar.Prod)
#GyrPressure[, 1] <- probOcc[, 1]

# Loop to fill out the rest of the years
for (t in 1:(N_years - 1)) {
  for (x in 1:N_areas) {
    
    # Current total density
    totalDensity[x, t] <- AdultDensity[x, t] + JuvenileDensity[x, t]
    
    # Predict next year's occupancy 
    probOcc[x, t + 1] <- plogis(alphaPtar.Occ[x] + betaPtar.Occ * totalDensity[x, t])
    
    # Calculate gyrpressure (average of current and next occupancy, or lagged only)
    GyrPressure[x, t] <- 0.5 * probOcc[x, t] + 0.5 * probOcc[x, t + 1]
    if (t == 1) {
      GyrPressure[x, t] <- probOcc[x, 1]  # use initial occupancy for first step
    } else {
      #GyrPressure[x, t] <- probOcc[x, t - 1]
      GyrPressure[x, t] <- 0.5 * probOcc[x, t-1] + 0.5 * probOcc[x, t]
    }
    
    # 4. Survival for next year (logistic, depends on g_index)
    S[x, t + 1] <- plogis(qlogis(Mu.S[x]) + betaGyr.S * GyrPressure[x, t])
    
    # 5. Recruitment for next year
    terrProd[x, t + 1] <- exp(alphaPtar.Prod[x] + betaPtar.Prod * totalDensity[x, t])
    # R_year[x, t + 1] <- exp(log(Mu.R[x]) + betaR.R * terrProd[x, t + 1])
    R_year[x, t + 1] <- exp(log(Mu.R[x])+ betaR.R * RodentOcc[x, t])
    
    # 6. Update densities for next year
    AdultDensity[x, t + 1] <- totalDensity[x, t] * S[x, t + 1]
    JuvenileDensity[x, t + 1] <- AdultDensity[x, t + 1] * R_year[x, t + 1]
    }
  }


matplot(t(probOcc), type='l', lty=1, main="Gyrfalcon Occupancy", ylab="Probability", xlab="Year")
matplot(t(AdultDensity), type='l', lty=1, main="Ptarmigan Adult Density", ylab="Density", xlab="Year")
matplot(t(JuvenileDensity), type='l', lty=1, main="Ptarmigan Juvenile Density", ylab="Density", xlab="Year")
matplot(t(totalDensity), type='l', lty=1, main="Ptarmigan total Density", ylab="Density", xlab="Year")
matplot(t(R_year), type='l', lty=1, main="Ptarmigan Recruitment", ylab="Recruitment", xlab="Year")
matplot(t(S), type='l', lty=1, main="Ptarmigan Survival", ylab="Survival", xlab="Year")
matplot(t(RodentOcc), type='l', lty=1, main="Rodent Occupancy", ylab="Rodent occupancy standardized", xlab="Year")




