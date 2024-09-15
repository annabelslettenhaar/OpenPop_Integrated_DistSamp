library(tidyverse)
library(coda)

# SETUP #
#-------#

## Set seed
mySeed <- 32
set.seed(mySeed)

## Set number of chains
nchains <- 5

## Set min and max years
minYear <- 2007
maxYear <- 2021

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')

## Create plotting directory
dir.create("Plots")

# RETRIEVE INPUT DATA AND POSTERIOR SAMPLES #
#-------------------------------------------#

## Read in input data
input_data <- readRDS("RypeData_forIM.rds")
LT_data <- readRDS("LT_data.rds")

## Reinstate model toggles
toggles <- readRDS("ModelRunToggles.rds")
list2env(toggles, globalenv())

## Retrieve seeds
seedList <- read.delim("inputSeeds.txt", header = FALSE)

## Assemble separate chain runs into one mcmc list
IDSM.out <- coda::mcmc.list()
for(i in 1:nchains){
  originSeed <- seedList[i, 1]
  chainSeed <- seedList[i, 2]
  IDSM.out[[i]] <- readRDS(paste0("rypeIDSM_dHN_multiArea_realData_allAreas_GNU_", originSeed, "_", chainSeed, ".rds"))
}

# TIDY UP POSTERIOR SAMPLES #
#---------------------------#

IDSM.out.tidy <- tidySamples(IDSM.out = IDSM.out,
                             save = TRUE,
                             fileName = "rypeIDSM_dHN_multiArea_realData_allAreas_tidy.rds")

#IDSM.out.tidy <- readRDS("rypeIDSM_dHN_multiArea_realData_allAreas_tidy.rds")

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

#PostSum.list <- readRDS("PosteriorSummaries_byArea.rds")

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

plotPosteriorDens_VR(mcmc.out = IDSM.out.tidy,
                     N_areas = input_data$nim.constant$N_areas,
                     area_names = input_data$nim.constant$area_names,
                     N_years = input_data$nim.constant$N_years,
                     minYear = minYear,
                     survAreaIdx = input_data$nim.constants$SurvAreaIdx,
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
                    meanCov = input_data$nim.constants$RodentOcc_meanCov,
                    sdCov = input_data$nim.constants$RodentOcc_sdCov,
                    covData = input_data$nim.data$RodentOcc,
                    N_areas = input_data$nim.constant$N_areas,
                    area_names = input_data$nim.constant$area_names,
                    fitRodentCov = fitRodentCov)
}


# OPTIONAL: PLOT DETECTION FUNCTIONS #
#------------------------------------#

plotDetectFunction(mcmc.out = IDSM.out.tidy,
                   maxDist = input_data$nim.constants$W,
                   N_areas = input_data$nim.constant$N_areas,
                   area_names = input_data$nim.constant$area_names)

  
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
                                        areas = listAreas(), areaAggregation = TRUE)

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

# OPTIONAL: MODEL COMPARISON #
#----------------------------#

# plotModelComparison(modelPaths = c("rypeIDSM_dHN_multiArea_realData_allAreas_tidy.rds",
#                                    "rypeIDSM_dHN_multiArea_realData_allAreas_noTelemetry_tidy.rds"),
#                     modelChars = c("Including telemetry",
#                                    "Without telemetry"),
#                     N_areas = input_data$nim.constants$N_areas,
#                     area_names = input_data$nim.constant$area_names,
#                     N_sites = input_data$nim.constants$N_sites,
#                     N_years = input_data$nim.constants$N_years,
#                     minYear = minYear,
#                     maxYear = maxYear,
#                     max_years = input_data$nim.constants$max_years,
#                     survAreaIdx = input_data$nim.constants$SurvAreaIdx,
#                     plotPath = "Plots/Comp_noTelemetry",
#                     returnData = FALSE)
