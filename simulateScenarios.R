library(tidyverse)
library(sf)
library(terra)
library(parallel)
library(nimble)
library(coda)


# Limits and constants #
#----------------------#

set.seed(123) 
  
# N_areas <- input_data$nim.constants$N_areas
# N_ageC <- input_data$nim.constants$N_ageC
# N_years <- input_data$nim.constants$N_years
N_areas <- 3
N_ageC <- 2
N_years <- 30

# L <- input_data$nim.data$L
# W <- input_data$nim.constants$W
# pi <- 3.141593
# A <- input_data$nim.data$A

# totDens_meanCov <- input_data$nim.constants$totDens_meanCov
# totDens_sdCov <- input_data$nim.constants$totDens_sdCov


# Rodent data #
#-------------#

# Second order auto regressive model (AR2) to simulate rodent dynamics 

# simulateRodentAR2 <- function(N_areas, N_years, alpha, phi1, phi2, sigma) {
#   mat <- matrix(NA, nrow = N_areas, ncol = N_years)
#   for (i in 1:N_areas) {
#     mat[i, 1:2] <- rnorm(2, 0, 1)
#     for (t in 3:N_years) {
#       mat[i, t] <- alpha + phi1 * mat[i, t-1] + phi2 * mat[i, t-2] + rnorm(1, 0, sigma)
#     }
#   }
#   mat
# }

# More detailed version

simulateRodentAR2 <- function(N_areas, N_years,
                                   alpha = 0, phi1 = 0.6, phi2 = -0.3, sigma = 0.2,
                                   peak_interval = 4, peak_jitter = 0, peak_size = 3.5,
                                   baseline = 0.1,
                                   mean_val = 5.39, sd_val = 9.36) {
  mat <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  for (i in 1:N_areas) {
    mat[i, 1:2] <- baseline + runif(2, 0, 0.2)  # start near baseline
    
    # Generate peak years with randomness
    peak_years <- seq(peak_interval, N_years, by = peak_interval) +
      sample(-peak_jitter:peak_jitter, length(seq(peak_interval, N_years, by = peak_interval)), replace = TRUE)
    peak_years <- peak_years[peak_years > 2 & peak_years <= N_years]  # keep valid years
    
    for (t in 3:N_years) {
      shock <- ifelse(t %in% peak_years, peak_size, 0)
      val <- alpha + phi1 * mat[i, t-1] + phi2 * mat[i, t-2] + rnorm(1, 0, sigma) + shock
      mat[i, t] <- max(val, baseline)  # enforce positivity
    }
  }
  
  # Transform back to original scale
  # mat_real <- mat * sd_val + mean_val
  return(mat)
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
chosen <- "weak"
p <- params[[chosen]]

# Run simulation
RodentOcc <- simulateRodentAR2(N_areas = 3, N_years = 100, 
                         alpha = p$alpha, phi1 = p$phi1, phi2 = p$phi2, sigma = p$sigma)


matplot(t(RodentOcc), type='l', lty=1, main=paste("Scenario:", chosen),
        ylab="Rodent index", xlab="Year")


## Get temperature data and resample

d_temp <- readRDS("data/weather/springtemp.rds")
d_temp <- d_temp$data

extra_years <- 70
rows <- nrow(d_temp)
cols <- ncol(d_temp)

future_data <- d_temp[, sample(1:cols, extra_years, replace = TRUE)]
SpringTemp <- cbind(d_temp, future_data)

matplot(t(SpringTemp), type = "l", lty = 1,
        main = "Simulated Spring Temperatures (100 years)",
        xlab = "Year", ylab = "Temperature")

# Helper function to get posterior medians #
# -----------------------------------------#

extractPostMedians <- function(modelOutput, paramNames) {
  samps <- as.matrix(modelOutput)
  median_list <- list()
  for (pname in paramNames) {
    param_cols <- grep(paste0("^", pname), colnames(samps), value = TRUE)
    medians <- apply(samps[, param_cols, drop = FALSE], 2, median)
    median_list[[pname]] <- medians
  }
  return(median_list)
}

# Load model output
model_output <- readRDS("/cloud/project/rypeIDSM_dHN_gyrData_11-11_longrun_fullloop_Rodent_Temp.rds")


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

betaR.R <- post_ptar$betaR.R


# Population density #
#--------------------#

## Extract posterior median
post_pop <- extractPostMedians(modelOutput = model_output,
                              paramNames = c("Mu.D1"))


## Initial densities / population sizes

# Intercept
Mu.D1 <- post_pop$Mu.D1


# Putting models together #
#-------------------------#

# Starting values from theoretical example
alphaPtar.Occ <- -5
betaPtar.Occ <- 0.2
Mu.S <- 0.5
Mu.R <- 1.5
betaGyr.S <- -2 # -5 in example
betaR.R <- 0.3
betaTemp.R <- 0.1

sim.years <- 100 # Increase simulation years to make oscillations visible
N_years <- sim.years

# Initialize matrices
AdultDensity <- JuvenileDensity <- totalDensity <- matrix(NA, nrow = N_areas, ncol = N_years)
S <- matrix(NA, nrow = N_areas, ncol = N_years)
R_year <- matrix(NA, nrow = N_areas, ncol = N_years)
probOcc <- terrProd <- matrix(NA, nrow = N_areas, ncol = N_years)
GyrPressure <- matrix(NA, nrow = N_areas, ncol = N_years)

# Initial values for the first year
AdultDensity[, 1] <- Mu.D1 * 1000000 # Convert to individuals per square km instead of per square meter
JuvenileDensity[, 1] <- AdultDensity[, 1]*Mu.R
probOcc[, 1] <- plogis(alphaPtar.Occ)

AdultDensity[, 1] <- 6
JuvenileDensity[, 1] <- 2
probOcc[, 1] <- 0.5

terrProd[, 1] <- exp(alphaPtar.Prod)
totalDensity[, 1] <- AdultDensity[, 1] + JuvenileDensity[, 1]


# Loop to fill out the rest of the years
for (t in 1:(N_years - 1)) {
  for (x in 1:N_areas) {
    
    # Current total density
    totalDensity[x, t] <- AdultDensity[x, t] + JuvenileDensity[x, t]
    
    # Predict next year's occupancy 
    probOcc[x, t + 1] <- plogis(alphaPtar.Occ[x] + betaPtar.Occ * totalDensity[x, t])
    
    # Calculate gyrpressure (average of current and next occupancy, or lagged only)
    # GyrPressure[x, t] <- 0.5 * probOcc[x, t] + 0.5 * probOcc[x, t + 1]
    if (t == 1) {
      GyrPressure[x, t] <- probOcc[x, 1]  # use initial occupancy for first step
    } else {
      # GyrPressure[x, t] <- probOcc[x, t - 1]
      GyrPressure[x, t] <- 0.5 * probOcc[x, t-1] + 0.5 * probOcc[x, t]
    }
    
    # 4. Survival for next year (logistic, depends on g_index)
    S[x, t + 1] <- plogis(qlogis(Mu.S[x]) + betaGyr.S * GyrPressure[x, t])
    
    # 5. Recruitment for next year
    # terrProd[x, t + 1] <- exp(alphaPtar.Prod[x] + betaPtar.Prod * totalDensity[x, t])
    # R_year[x, t + 1] <- exp(log(Mu.R[x]))
    R_year[x, t + 1] <- exp(log(Mu.R[x]) + betaR.R * RodentOcc[x, t])
    # R_year[x, t + 1] <- exp(log(Mu.R[x]) + betaTemp.R * SpringTemp[x, t])
    # R_year[x, t + 1] <- exp(log(Mu.R[x]) + betaR.R * RodentOcc[x, t] + betaTemp.R * SpringTemp[x, t])
    
    # 6. Update densities for next year
    AdultDensity[x, t + 1] <- totalDensity[x, t] * S[x, t + 1]
    JuvenileDensity[x, t + 1] <- AdultDensity[x, t + 1] * R_year[x, t + 1]
    }
  }

par(mfrow = c(1, 2))
matplot(t(probOcc), type='l', lty=1, main="Gyrfalcon Occupancy", ylab="Probability", xlab="Year")
matplot(t(totalDensity), type='l', lty=1, main="Ptarmigan total Density", ylab="Density", xlab="Year")


# Two panels side by side
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# Gyrfalcon Occupancy
matplot(t(probOcc), type = 'l', lty = 1, lwd = 3, col = c("steelblue", "darkgreen", "firebrick"),
        main = "Gyrfalcon brood initiation", ylab = "Probability", xlab = "Time step")
#legend("bottomright", legend = c("Area 1", "Area 2", "Area 3"),
#       col = c("steelblue", "darkgreen", "firebrick"), lty = 1, lwd = 3, cex = 0.9)

# Ptarmigan Density
matplot(t(totalDensity), type = 'l', lty = 1, lwd = 3, col = c("steelblue", "darkgreen", "firebrick"),
        main = "Ptarmigan total density", ylab = "Density (individuals per sqkm)", xlab = "Time step")
#legend("bottomright", legend = c("Area 1", "Area 2", "Area 3"),
#       col = c("steelblue", "darkgreen", "firebrick"), lty = 1, lwd = 3, cex = 0.9)




matplot(t(AdultDensity), type='l', lty=1, main="Ptarmigan Adult Density", ylab="Density", xlab="Year")
matplot(t(JuvenileDensity), type='l', lty=1, main="Ptarmigan Juvenile Density", ylab="Density", xlab="Year")
matplot(t(R_year), type='l', lty=1, main="Ptarmigan Recruitment", ylab="Recruitment", xlab="Year")
matplot(t(S), type='l', lty=1, main="Ptarmigan Survival", ylab="Survival", xlab="Year")
matplot(t(RodentOcc), type='l', lty=1, main="Rodent Occupancy", ylab="Rodent occupancy standardized", xlab="Year")




