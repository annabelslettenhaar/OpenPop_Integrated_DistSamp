## Load packages required to define the pipeline
library(targets)
library(nimble)
library(nimbleDistance)
library(tidyverse)
library(ggplot2)
library(LivingNorwayR)
library(terra)
library(sf)
library(tmap)
library(parallel)

## Set seed
mySeed <- 32
set.seed(mySeed)

nchains <- 5

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')

## Set localities and time period of interest
#localities <- unname(read.csv("localities.csv", encoding = "utf-8"))[,1]
minYear <- 2007
maxYear <- 2021


## Set model toggles
areaAggregation <- TRUE # Area- vs. locality aggregation
R_perF <- FALSE # Recruitment per adult or per adult female
R_parent_drop0 <- TRUE # Drop observations of juveniles with no adults present
sumR.Level <- "line" # Aggregation level for reproduction data (line vs. group)
survVarT <- TRUE # Time variation in survival
fitRodentCov <- TRUE # Rodent covariate on reproduction
telemetryData <- TRUE # Include telemetry data

parallelMCMC <- TRUE

## Set target-specific options such as packages.
tar_option_set(packages = c("LivingNorwayR", "tidyverse", "qs"),
               format = "qs",
               memory = "transient", 
               garbage_collection = TRUE)

## Define Targets List
list(
  tar_target(
    Rype_arkiv,
    downloadLN(datasets = c("Fjellstyrene", "Statskog", "FeFo"), 
               versions = c(1.7, 1.8, 1.12), 
               save = TRUE)
  ),
  
  tar_target(
    duplTransects,
    listDuplTransects()
  ),
  
  tar_target(
    localities,
    listLocations()
  ),
  
  tar_target(
    areas,
    listAreas()
  ),
  
  tar_target(
    LT_data,
    wrangleData_LineTrans(DwC_archive_list = Rype_arkiv,
                          duplTransects = duplTransects,
                          #localities = localities,
                          areas = areas,
                          areaAggregation = areaAggregation,
                          minYear = minYear, maxYear = maxYear)
  ),
  
  tar_target(
    d_cmr,
    wrangleData_CMR(minYear = minYear)
  ),
  
  tar_target(
    d_rodent,
    wrangleData_Rodent(duplTransects = duplTransects,
                       #localities = localities,
                       areas = areas,
                       areaAggregation = TRUE,
                       minYear = minYear, maxYear = maxYear)
  ),
  
  tar_target(
    input_data,
    prepareInputData(d_trans = LT_data$d_trans, 
                     d_obs = LT_data$d_obs,
                     d_cmr = d_cmr,
                     d_rodent = d_rodent,
                     #localities = localities,
                     areas = areas,
                     areaAggregation = areaAggregation,
                     excl_neverObs = TRUE,
                     R_perF = R_perF,
                     R_parent_drop0 = R_parent_drop0,
                     sumR.Level = sumR.Level,
                     dataVSconstants = TRUE,
                     save = TRUE)
  ),
  
  tar_target(
    modelCode,
    writeModelCode(survVarT = survVarT,
                   telemetryData = telemetryData)
  ),
  
  
  tar_target(
    MCMC.seeds,
    expandSeed_MCMC(seed = mySeed, 
                    nchains = nchains)
  ),
  
  tar_target(
    model_setup,
    setupModel(modelCode = modelCode,
               R_perF = R_perF,
               nim.data = input_data$nim.data,
               nim.constants = input_data$nim.constants,
               survVarT = survVarT,
               fitRodentCov = fitRodentCov,
               testRun = TRUE, 
               nchains = nchains,
               initVals.seed = MCMC.seeds)
  ),
  
  tar_target(
    IDSM.out,
    runModel(parallelMCMC = parallelMCMC,
             model_setup = model_setup, 
             MCMC.seeds = MCMC.seeds, 
             input_data = input_data)
  ),
  
  tar_target(
    IDSM.out.tidy,
    tidySamples(IDSM.out = IDSM.out,
                save = TRUE,
                fileName = "rypeIDSM_dHN_multiArea_realData_allAreas_tidy.rds")
  ),
  
  tar_target(
    PostSum.list,
    summarisePost_areas(mcmc.out = IDSM.out.tidy, 
                        N_areas = input_data$nim.constant$N_areas, 
                        area_names = input_data$nim.constant$area_names, 
                        N_sites = input_data$nim.constant$N_sites, 
                        min_years = input_data$nim.constant$min_years, 
                        max_years = input_data$nim.constant$max_years, 
                        minYear = minYear, maxYear = maxYear,
                        fitRodentCov = fitRodentCov,
                        save = TRUE)
  ),
  
  tar_target(
    mcmc.tracePlots,
    plotMCMCTraces(mcmc.out = IDSM.out.tidy,
                   fitRodentCov = fitRodentCov,
                   survVarT = survVarT),
    format = "file"
  ),
  
  tar_target(
    time.seriesPlots,
    plotTimeSeries(mcmc.out = IDSM.out.tidy, 
                   N_areas = input_data$nim.constant$N_areas, 
                   area_names = input_data$nim.constant$area_names, 
                   N_sites = input_data$nim.constant$N_sites, 
                   min_years = input_data$nim.constant$min_years, 
                   max_years = input_data$nim.constant$max_years, 
                   minYear = minYear, maxYear = maxYear,
                   VitalRates = TRUE, DetectParams = TRUE, Densities = TRUE),
    format = "file"
  ),

  tar_target(
    post.densPlots,
    plotPosteriorDens_VR(mcmc.out = IDSM.out.tidy,
                         N_areas = input_data$nim.constant$N_areas, 
                         area_names = input_data$nim.constant$area_names, 
                         N_years = input_data$nim.constant$N_years,
                         minYear = minYear,
                         survAreaIdx = input_data$nim.constants$SurvAreaIdx,
                         survVarT = survVarT,
                         fitRodentCov = fitRodentCov) 
  ),
  
  tar_target(
    cov.predPlots,
    plotCovPrediction(mcmc.out = IDSM.out.tidy,
                      effectParam = "betaR.R",
                      covName = "Rodent occupancy",
                      minCov = 0, 
                      maxCov = 1,
                      meanCov = d_rodent$meanCov,
                      sdCov = d_rodent$sdCov,
                      N_areas = input_data$nim.constant$N_areas, 
                      area_names = input_data$nim.constant$area_names,
                      fitRodentCov = fitRodentCov)
  ),
  
  tar_target(
    DD.checks,
    checkDD(mcmc.out = IDSM.out.tidy, 
            N_areas = input_data$nim.constant$N_areas, 
            area_names = input_data$nim.constant$area_names, 
            N_sites = input_data$nim.constant$N_sites, 
            min_years = input_data$nim.constant$min_years, 
            max_years = input_data$nim.constant$max_years),
    format = "file"
  ),
  
  tar_target(
    VRcorr.checks,
    checkVRcorrs(mcmc.out = IDSM.out.tidy, 
                 N_areas = input_data$nim.constant$N_areas, 
                 area_names = input_data$nim.constant$area_names, 
                 area_coord = LT_data$d_coord,
                 min_years = input_data$nim.constant$min_years, 
                 max_years = input_data$nim.constant$max_years)
  ),
  
  tar_target(
    VarDecomp,
    plotVarDecomposition(mcmc.out = IDSM.out.tidy, 
                         N_areas = input_data$nim.constants$N_areas, 
                         N_years = input_data$nim.constants$N_years, 
                         fitRodentCov = fitRodentCov, 
                         RodentOcc_data = input_data$nim.data$RodentOcc)
  ),
  
  tar_target(
    NorwayMunic.map,
    setupMap_NorwayMunic(shp.path = "data/norway_municipalities/norway_municipalities.shp",
                         d_trans = LT_data$d_trans,
                         areas = areas, areaAggregation = areaAggregation)
  ),
  
  tar_target(
    mapPlots,
    plotMaps(PostSum.list = PostSum.list, 
             mapNM = NorwayMunic.map,
             minYear = minYear, maxYear = maxYear,
             fitRodentCov = fitRodentCov),
    format = "file"
  ),
  
  tar_target(
    latitudePlots,
    plotLatitude(PostSum.list = PostSum.list, 
                 area_coord = LT_data$d_coord,
                 minYear = minYear, maxYear = maxYear,
                 fitRodentCov = fitRodentCov),
    format = "file"
  ),
  
  tar_target(
    GT_estimates,
    extract_GenTime(mcmc.out = IDSM.out.tidy, 
                    N_areas = input_data$nim.constants$N_areas, 
                    area_names = input_data$nim.constant$area_names, 
                    area_coord = LT_data$d_coord,
                    mapNM = NorwayMunic.map,
                    save = TRUE)
  )
)

# Test by running tar_manifest(fields = all_of("command")) and tar_visnetwork() in the console

# Run workflow using tar_make() in the console
# Check with tar_network() in the console

# We can then use tar_read() and tar_load() to inspect and work with results
