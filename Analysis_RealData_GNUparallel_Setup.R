library(tidyverse)
library(nimble)

# SETUP #
#-------#

## Set seed
mySeed <- 32
set.seed(mySeed)

## Set number of chains
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


## Set and store switches/toggles 

# (Re-)downloading data
# downloadData <- FALSE
downloadData <- TRUE

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
telemetryData <- TRUE

# Test run or not
testRun <- FALSE

# Run MCMC in parallel
parallelMCMC <- FALSE

# Store as .RDS
storeToggles(downloadData = downloadData,
             R_perF = R_perF,
             R_parent_drop0 = R_parent_drop0,
             sumR.Level = sumR.Level,
             survVarT = survVarT,
             fitRodentCov = fitRodentCov,
             telemetryData = telemetryData,
             testRun = testRun,
             parallelMCMC = parallelMCMC)


# DOWNLOAD/FETCH DATA #
#---------------------#

if(downloadData){
  #Rype_arkiv <- downloadLN(datasets = "Fjellstyrene", versions = 1.6, save = TRUE)
  Rype_arkiv <- downloadLN(datasets = c("Fjellstyrene", "Statskog", "FeFo"), versions = c(1.7, 1.8, 1.12), save = TRUE)
}else{
  stop("downloadData = FALSE not supported yet. There is an issue with encoding when using LivingNorwayR::initializeDwCArchive() that needs to be resolved first.")
  #Rype_arkiv <- initializeDwCArchive("data/Rype_arkiv.zip")
}


# WRANGLE LINE TRANSECT DATA #
#----------------------------#

## Set localities/areas and time period of interest
localities <- listLocations()
areas <- listAreas()
#areas <- listAreas()[c(5, 17, 34)]
minYear <- 2007
maxYear <- 2021

## List duplicate transects to remove
duplTransects <- listDuplTransects()

## Extract transect and observational data from DwC archive
LT_data <- wrangleData_LineTrans(DwC_archive_list = Rype_arkiv, 
                                 duplTransects = duplTransects,
                                 #localities = localities,
                                 areas = areas,
                                 areaAggregation = TRUE,
                                 minYear = minYear, maxYear = maxYear)

## Save line transect data for use in post-processing
saveRDS(LT_data, file = "LT_data.rds")

# WRANGLE KNOWN FATE CMR DATA #
#-----------------------------#

## Read in and reformat CMR data
d_cmr <- wrangleData_CMR(minYear = minYear)


# WRANGLE RODENT DATA #
#---------------------#

## Load and reformat rodent data
d_rodent <- wrangleData_Rodent(duplTransects = duplTransects,
                               #localities = localities,
                               areas = areas,
                               areaAggregation = TRUE,
                               minYear = minYear, maxYear = maxYear)


# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Reformat data into vector/array list for analysis with Nimble
input_data <- prepareInputData(d_trans = LT_data$d_trans, 
                               d_obs = LT_data$d_obs,
                               d_cmr = d_cmr,
                               d_rodent = d_rodent,
                               #localities = localities, 
                               areas = areas,
                               areaAggregation = TRUE,
                               excl_neverObs = TRUE,
                               R_perF = R_perF,
                               R_parent_drop0 = R_parent_drop0,
                               sumR.Level = "line",
                               dataVSconstants = TRUE,
                               save = TRUE)


# SEED SIMULATION #
#-----------------#

## Expand seeds for simulating initial values
MCMC.seeds <- expandSeed_MCMC(seed = mySeed, 
                              nchains = nchains)

# --> write seed combinations into a .txt file (tab-delimited)

seedList <- list()
seedList <- c()

for(i in 1:length(MCMC.seeds)){
  #seedList[[i]] <- paste0(mySeed, "\t", MCMC.seeds[i]) 
  seedList <- c(seedList, paste0(mySeed, "\t", MCMC.seeds[i]) )
}

writeLines(seedList, 
           "inputSeeds.txt")
