library(tidyverse)
library(sf)
library(terra)
library(parallel)
library(nimble)
library(coda)


# SETUP #
#-------#

## Set seed
mySeed <- 83
set.seed(mySeed)

## Set number of chains, iterations, burn in and thinning
nchains <- 3
niter <- 200000
nburn <- 140000
nthin <- 20

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')


## Set and store switches/toggles 

# Aggregation to area level
areaAggregation <- TRUE 

# Recruitment per adult or per adult female
R_perF <- FALSE

# Drop observations of juveniles with no adults present
R_parent_drop0 <- TRUE

# Aggregation level for reproduction data
# NOTE: if this is not defined, will default to group level
sumR.Level <- "line" # Summing at the line level

# Time variation in survival
survVarT <- TRUE

# Rodent covariate on reproduction
fitRodentCov <- TRUE

# Use of telemetry data from Lierne
telemetryData <- FALSE

# Test run or not
testRun <- TRUE

# Run MCMC in parallel
parallelMCMC <- FALSE

# Fully closed loop in predator-prey model
fullLoopPP <- TRUE


# WRANGLE LINE TRANSECT DATA #
#----------------------------#

## Set localities/areas and time period of interest
localities <- listLocations()
areas <- c("Hardangervidda", 
           "Dovrefjell", 
           "Børgefjell")
minYear <- 1991
maxYear <- 2020

## List duplicate transects to remove
duplTransects <- listDuplTransects()

## Extract transect and observational data from DwC archive
LT_data <- wrangleData_DwCPtar(#localities = localities,
  areas = areas,
  areaAggregation = areaAggregation,
  minYear = minYear, maxYear = maxYear)



# WRANGLE RODENT DATA #
#---------------------#

## Load and reformat rodent data
d_rodent <- wrangleData_RodentGyr(#localities = localities,
                                  areas = areas,
                                  areaAggregation = areaAggregation,
                                  minYear = minYear, maxYear = maxYear)

# WRANGLE WEATHER DATA #
#----------------------#

sourceDir('R_weather')

# Using gyrfalcon territory coordinates
weather_data <- wrangleData_Weather_GT(areas = areas,
                                       areaAggregation = areaAggregation,
                                       minYear = minYear,
                                       maxYear = maxYear)

# The weather_data above are not used in the actual analysis: can be cleaned up later


## Using ptarmigan line transect coordinates

# d_temp <- wrangleData_Temp(minYear = minYear,
#                            maxYear = maxYear,
#                            areas = areas,
#                            startday = 121,
#                            endday = 153)
# saveRDS(d_temp, "data/weather/springtemp.rds")

d_temp <- readRDS("data/weather/springtemp.rds")



# WRANGLE GYRFALCON DATA #
#------------------------#

## Load gyr pressure data
d_gyr <- wrangleData_GyrPressure(#localities = localities,
  areas = areas,
  areaAggregation = areaAggregation,
  minYear = minYear, maxYear = maxYear)

## Load gyr productivity data
d_gyrprod <- wrangleData_GyrProd_nested(minYear, 
                                        maxYear)

## Load gyr occupancy data
d_gyrocc <- wrangleData_GyrOcc_agg(minYear,
                                   maxYear)

# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Define mean and sd for standardizing ptarmigan density in the model (per area)
totDens_meanCov <- c(1.4e-05, 2.1e-05, 2.5e-05)
totDens_sdCov <- c(8e-06, 8e-06, 1.7e-05)

# ## Define mean and sd for standardizing GyrPressure in the model (occ + prod, indiv numbers)
# GyrPressure_meanCov <- c(18.5, 12.3, 12.1)
# GyrPressure_sdCov <- c(6.43, 3.91, 3.43)

## Define mean and sd for standardizing GyrPressure in the model (only occ probability)
GyrPressure_meanCov <- c(0.371, 0.35, 0.283)
GyrPressure_sdCov <- c(0.1136, 0.0742, 0.0716)

## Reformat data into vector/array list for analysis with Nimble
input_data <- prepareInputData_Integ(d_trans = LT_data$d_trans, 
                                     d_obs = LT_data$d_obs,
                                     #d_cmr = d_cmr,
                                     d_rodent = d_rodent,
                                     d_gyr = d_gyr, # GyrPressure covariate
                                     d_gyrocc = d_gyrocc, # Occupancy data
                                     d_gyrprod = d_gyrprod, # Productivity data
                                     d_SD = weather_data$d_SD, # For analysis including weather
                                     d_temp = d_temp,
                                     #localities = localities, 
                                     areas = areas,
                                     areaAggregation = areaAggregation,
                                     excl_neverObs = TRUE,
                                     R_perF = R_perF,
                                     R_parent_drop0 = R_parent_drop0,
                                     sumR.Level = "line",
                                     totDens_meanCov = totDens_meanCov,
                                     totDens_sdCov = totDens_sdCov, 
                                     GyrPressure_meanCov = GyrPressure_meanCov,
                                     GyrPressure_sdCov = GyrPressure_sdCov,
                                     fullLoopPP = fullLoopPP,
                                     dataVSconstants = TRUE,
                                     save = TRUE)


# MODEL SETUP #
#-------------#

## Write model code
modelCode <- writeModelCode_Integ(survVarT = survVarT,
                                  telemetryData = telemetryData,
                                  fullLoopPP = fullLoopPP)

## Expand seeds for simulating initial values
MCMC.seeds <- expandSeed_MCMC(seed = mySeed, 
                              nchains = nchains)

#MCMC.seeds <- MCMC.seeds[1]

## Setup for model using nimbleDistance::dHN
model_setup <- setupModel_Integ(modelCode = modelCode,
                                R_perF = R_perF,
                                survVarT = survVarT, 
                                fitRodentCov = fitRodentCov,
                                fullLoopPP = fullLoopPP,
                                nim.data = input_data$nim.data,
                                nim.constants = input_data$nim.constants,
                                testRun = testRun, 
                                nchains = nchains,
                                #nchains = 1,
                                niter = niter,
                                nburn = nburn,
                                nthin = nthin,
                                initVals.seed = MCMC.seeds)


# MODEL (TEST) RUN #
#------------------#

if(!parallelMCMC){
  t.start <- Sys.time()
  IDSM.out <- nimbleMCMC(code = model_setup$modelCode,
                         data = input_data$nim.data, 
                         constants = input_data$nim.constants,
                         inits = model_setup$initVals, 
                         monitors = model_setup$modelParams,
                         nchains = model_setup$mcmcParams$nchains, 
                         niter = model_setup$mcmcParams$niter, 
                         nburnin = model_setup$mcmcParams$nburn, 
                         thin = model_setup$mcmcParams$nthin, 
                         samplesAsCodaMCMC = TRUE, 
                         setSeed = MCMC.seeds)
  Sys.time() - t.start
  
  
}else{
  
  ## Add toggles to constants
  input_data$nim.constants$fitRodentCov <- fitRodentCov
  input_data$nim.constants$survVarT <- survVarT
  input_data$nim.constants$R_perF <- R_perF
  input_data$nim.constants$telemetryData <- telemetryData
  
  ## Set up cluster
  this_cluster <- makeCluster(model_setup$mcmcParams$nchains)
  #clusterEvalQ(this_cluster, library(nimble))
  #clusterEvalQ(this_cluster, library(nimbleDistance))
  
  ## Collect chain-specific information
  per_chain_info <- vector("list", model_setup$mcmcParams$nchains)
  for(i in 1:model_setup$mcmcParams$nchains){
    per_chain_info[[i]] <- list(mySeed = MCMC.seeds[i],
                                inits = model_setup$initVals[[i]])
  }
  
  ## Run chains in parallel
  t.start <- Sys.time()
  IDSM.out <- parLapply(cl = this_cluster, 
                        X = per_chain_info, 
                        fun = runMCMC_allcode, 
                        model_setup = model_setup,
                        input_data = input_data)
  Sys.time() - t.start
  
  
  stopCluster(this_cluster)
  
}

saveRDS(IDSM.out, file = "rypeIDSM_dHN_gyrData_11-11_longrun_fullloop_Rodent_Temp.rds")


# # TIDY UP POSTERIOR SAMPLES #
# #---------------------------#
# 
# IDSM.out.tidy <- tidySamples(IDSM.out = IDSM.out, 
#                              save = FALSE,
#                              fileName = "rypeIDSM_dHN_multiArea_gyrData_gyrCov3_tidy.rds")



# MAKE POSTERIOR SUMMARIES PER AREA #
#-----------------------------------#

# PostSum.list <- summarisePost_areas_gyr_integ(mcmc.out = onlyOcc, 
#                                               N_areas = input_data$nim.constant$N_areas, 
#                                               area_names = input_data$nim.constant$area_names, 
#                                               N_sites = input_data$nim.constant$N_sites, 
#                                               min_years = input_data$nim.constant$min_years, 
#                                               max_years = input_data$nim.constant$max_years, 
#                                               minYear = minYear, maxYear = maxYear,
#                                               fitRodentCov = fitRodentCov,
#                                               save = FALSE)
# Above code can be cleaned up later

# Adjusted version

var_list <- c("totDens_raw", "probOcc", "terrProd",
              "GyrPressure_raw", "R_year", "S",
              
              "Mu.S", "Mu.R", "alphaPtar.Prod", "alphaPtar.Occ",
              
              "betaR.R", "betaPtar.Prod", "betaPtar.Occ", "betaGyr.S", "betaTemp.R")

var_type <- c("area_year", "area_year", "area_year",
              "area_year", "area_year", "area_year",
              
              "area", "area", "area", "area",
              
              "overall", "overall", "overall", "overall", "overall")

postSum_list <- postSum(samps = as.matrix(IDSM.out), 
                        var_list = var_list,
                        var_type = var_type,
                        N_areas = 3)

## Summarize results in a table

postSumTable <- createTablePosteriorSummaries(postSum_list)


# PLOT COVARIATE EFFECT SIZES #
#-----------------------------#

CovPlot <- visualiseBetaPosteriors(IDSM_out = IDSM.out,
                          cov_names = c("betaR.R", "betaPtar.Occ",
                                        "betaPtar.Prod", "betaGyr.S",
                                        "betaTemp.R"))

# PLOT ESTIMATE TIME SERIES #
#---------------------------#

TimeSeriesPlot <- plotTimeSeries_gyr(IDSM_out = IDSM.out,
                                     minYear = minYear)



















