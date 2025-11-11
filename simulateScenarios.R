
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

# ## Area-specific detection parameters
# h.mu.dd <- runif(1, 3.5, 5.5)
# h.sigma.dd <- runif(1, 0.05, 0.2)
# 
# #mu.dd <- rnorm(N_areas, h.mu.dd, sd = h.sigma.dd)
# mu.dd <- rep(h.mu.dd, N_areas)
# 
# # sigmaT.dd <- runif(1, 0.05, 0.2)
# sigmaR.dd <- runif(1, 0.05, 0.2)
# 
# sigma <- esw <- p <- matrix(NA, nrow = N_areas, ncol = N_years)
# 
# 
# # epsT.dd <- rep(0, N_years)
# #epsT.dd <- rnorm(N_years, 0, sd = sigmaT.dd)
# epsR.dd <- matrix(0, nrow = N_areas, ncol = N_years)
# #epsR.dd <- matrix(rnorm(N_areas*N_years, 0, sigmaR.dd), nrow = N_areas, ncol = N_years)
# 
# for(x in 1:N_areas){
#   # sigma[x, 1:N_years] <- exp(mu.dd[x] + epsT.dd[1:N_years] + epsR.dd[x, 1:N_years])
#   sigma[x, 1:N_years] <- exp(mu.dd[x] + epsR.dd[x, 1:N_years])
# }
# 
# for(x in 1:N_areas){
#   esw[x, 1:N_years] <- sqrt(pi * sigma[x, 1:N_years]^2 / 2) 
#   p[x, 1:N_years] <- min(esw[x, 1:N_years], W) / W
# }

# Something like this:
ptarDens <- posterior_samples$ptardens

# Or with noise
ptarDens <- ptarDens * exp(rnorm(length(ptarDens), 0, 0.1))


# Population model #
#------------------#

## Initial densities / population sizes
Mu.D1 <- rep(NA, N_areas) # Replace with model estimates
#sigma.D <- runif(N_areas, 0.1, 2)

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
  totDens_raw = totDens_raw,
  totDens_std = totDens_std,
  N_exp = N_exp,
  N_tot_exp = N_tot_exp,
  
  betaPtar.Occ = betaPtar.Occ,
  betaPtar.Prod = betaPtar.Prod,
  alphaPtar.Occ = alphaPtar.Occ,
  alphaPtar.Prod = alphaPtar.Prod,
  
  terrProd = terrProd, 
  probOcc = probOcc, 
  
  betaTemp.R = betaTemp.R,
  
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

