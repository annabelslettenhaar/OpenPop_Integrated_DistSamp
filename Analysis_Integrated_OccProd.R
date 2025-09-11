library(tidyverse)
library(sf)
library(terra)
library(parallel)
library(nimble)

# SETUP #
#-------#

## Set seed
mySeed <- 83
set.seed(mySeed)

## Set number of chains, iterations, burn in and thinning
nchains <- 3
niter <- 10000
nburn <- 6000
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
survVarT <- FALSE

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
weather_data <- wrangleData_Weather(areas = areas,
                                    areaAggregation = areaAggregation,
                                    minYear = minYear, maxYear = maxYear)


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

## Reformat data into vector/array list for analysis with Nimble
input_data <- prepareInputData_Integ(d_trans = LT_data$d_trans, 
                                     d_obs = LT_data$d_obs,
                                    #d_cmr = d_cmr,
                                     d_rodent = d_rodent,
                                     d_gyr = d_gyr, # GyrPressure covariate
                                     d_gyrocc = d_gyrocc, # Occupancy data
                                     d_gyrprod = d_gyrprod, # Productivity data
                                     d_SD = weather_data$d_SD, # For analysis including weather
                                    #localities = localities, 
                                     areas = areas,
                                     areaAggregation = areaAggregation,
                                     excl_neverObs = TRUE,
                                     R_perF = R_perF,
                                     R_parent_drop0 = R_parent_drop0,
                                     sumR.Level = "line",
                                     totDens_meanCov = totDens_meanCov,
                                     totDens_sdCov = totDens_sdCov, 
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

MCMC.seeds <- MCMC.seeds[1]

## Setup for model using nimbleDistance::dHN
model_setup <- setupModel_Integ(modelCode = modelCode,
                                R_perF = R_perF,
                                survVarT = survVarT, 
                                fitRodentCov = fitRodentCov,
                                fullLoopPP = fullLoopPP,
                                nim.data = input_data$nim.data,
                                nim.constants = input_data$nim.constants,
                                testRun = testRun, 
                                #nchains = nchains,
                                nchains = 1,
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

saveRDS(IDSM.out, file = "rypeIDSM_dHN_gyrData_occprod_08-09_allcov_oneslope.rds")


# TIDY UP POSTERIOR SAMPLES #
#---------------------------#

IDSM.out.tidy <- tidySamples(IDSM.out = IDSM.out, 
                             save = FALSE,
                             fileName = "rypeIDSM_dHN_multiArea_gyrData_gyrCov3_tidy.rds")



# MAKE POSTERIOR SUMMARIES PER AREA #
#-----------------------------------#

PostSum.list <- summarisePost_areas(mcmc.out = IDSM.out.tidy, 
                                    N_areas = input_data$nim.constant$N_areas, 
                                    area_names = input_data$nim.constant$area_names, 
                                    N_sites = input_data$nim.constant$N_sites, 
                                    min_years = input_data$nim.constant$min_years, 
                                    max_years = input_data$nim.constant$max_years, 
                                    minYear = minYear, maxYear = maxYear,
                                    fitRodentCov = fitRodentCov,
                                    save = TRUE)


# OPTIONAL: MCMC TRACE PLOTS #
#----------------------------#

plotMCMCTraces(mcmc.out = IDSM.out.tidy,
               fitRodentCov = fitRodentCov,
               survVarT = survVarT)


# OPTIONAL: TIME SERIES PLOTS #
#-----------------------------#

plotTimeSeries(mcmc.out = IDSM.out.tidy, 
               N_areas = input_data$nim.constant$N_areas, 
               area_names = input_data$nim.constant$area_names, 
               N_sites = input_data$nim.constant$N_sites, 
               min_years = input_data$nim.constant$min_years, 
               max_years = input_data$nim.constant$max_years, 
               minYear = minYear, maxYear = maxYear,
               VitalRates = TRUE, DetectParams = TRUE, Densities = TRUE)


# OPTIONAL: PLOT VITAL RATE POSTERIORS #
#--------------------------------------#
# Needs to be adjusted to work without providing telemetry data
plotPosteriorDens_VR_Gyr(mcmc.out = IDSM.out.tidy,
                         N_areas = input_data$nim.constant$N_areas, 
                         area_names = input_data$nim.constant$area_names, 
                         N_years = input_data$nim.constant$N_years,
                         minYear = minYear,
                         #survAreaIdx = input_data$nim.constants$SurvAreaIdx,
                         survVarT = survVarT,
                         fitRodentCov = fitRodentCov) 


# OPTIONAL: PLOT COVARIATE PREDICTIONS #
#--------------------------------------#

if(fitRodentCov){
  plotCovPrediction(mcmc.out = IDSM.out.tidy,
                    effectParam = "betaR.R",
                    covName = "Rodent occupancy",
                    minCov = 0, 
                    maxCov = 1,
                    meanCov = d_rodent$meanCov,
                    sdCov = d_rodent$sdCov,
                    covData = d_rodent$rodentAvg,
                    N_areas = input_data$nim.constant$N_areas, 
                    area_names = input_data$nim.constant$area_names,
                    fitRodentCov = fitRodentCov)
}


# OPTIONAL: CHECK WITHIN-AREA DENSITY DEPENDENCE #
#------------------------------------------------#

checkDD(mcmc.out = IDSM.out.tidy, 
        N_areas = input_data$nim.constant$N_areas, 
        area_names = input_data$nim.constant$area_names, 
        N_sites = input_data$nim.constant$N_sites, 
        min_years = input_data$nim.constant$min_years, 
        max_years = input_data$nim.constant$max_years)


# OPTIONAL: CHECK VITAL RATE SAMPLING CORRELATIONS #
#--------------------------------------------------#

checkVRcorrs(mcmc.out = IDSM.out.tidy, 
             N_areas = input_data$nim.constant$N_areas, 
             area_names = input_data$nim.constant$area_names, 
             area_coord = LT_data$d_coord,
             min_years = input_data$nim.constant$min_years, 
             max_years = input_data$nim.constant$max_years)


# OPTIONAL: CALCULATE AND PLOT VARIANCE DECOMPOSITION #
#-----------------------------------------------------#

plotVarDecomposition(mcmc.out = IDSM.out.tidy, 
                     N_areas = input_data$nim.constants$N_areas, 
                     N_years = input_data$nim.constants$N_years, 
                     fitRodentCov = fitRodentCov, 
                     RodentOcc_data = input_data$nim.data$RodentOcc,
                     saveResults = TRUE)


# OPTIONAL: MAP PLOTS #
#---------------------#

## Make map of Norwegian municipalities ("fylke")
NorwayMunic.map <- setupMap_NorwayMunic(shp.path = "data/norway_municipalities/norway_municipalities.shp",
                                        d_trans = LT_data$d_trans,
                                        areas = areas, areaAggregation = areaAggregation)

## Plot population growth rate, density, and vital rates on map
plotMaps(PostSum.list = PostSum.list, 
         mapNM = NorwayMunic.map,
         minYear = minYear, maxYear = maxYear,
         fitRodentCov = fitRodentCov)


# OPTIONAL: LATITUDE PATTERN PLOTS #
#----------------------------------#

plotLatitude(PostSum.list = PostSum.list, 
             area_coord = LT_data$d_coord,
             minYear = minYear, maxYear = maxYear,
             fitRodentCov = fitRodentCov)


# OPTIONAL: GENERATION TIME #
#---------------------------#

GT_estimates <- extract_GenTime(mcmc.out = IDSM.out.tidy, 
                                N_areas = input_data$nim.constants$N_areas, 
                                area_names = input_data$nim.constant$area_names, 
                                area_coord = LT_data$d_coord,
                                mapNM = NorwayMunic.map,
                                save = TRUE)


# OPTIONAL: MODEL COMPARISON (PLOTS) #
#------------------------------------#

plotModelComparison(modelPaths = c("rypeIDSM_dHN_multiArea_realData_allAreas_tidy.rds",
                                   "rypeIDSM_dHN_multiArea_realData_allAreas_tidy_noTelemetry.rds"), 
                    modelChars = c("Including telemetry",
                                   "Without telemetry"), 
                    N_areas = input_data$nim.constants$N_areas, 
                    area_names = areas, 
                    N_sites = input_data$nim.constants$N_sites, 
                    N_years = input_data$nim.constants$N_years, 
                    minYear = minYear, 
                    maxYear = maxYear, 
                    max_years = input_data$nim.constants$max_years, 
                    survAreaIdx = input_data$nim.constants$SurvAreaIdx, 
                    plotPath = "Plots/Comp_noTelemetry", 
                    returnData = FALSE)
