library(tidyverse)

# SETUP #
#-------#

## Source all functions in "R" folder
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir('R')


## Set switches 

# (Re-)download data (TRUE) or use existing data in 'data' folder (FALSE)
downloadData <- TRUE

# NI: estimate parameters by pre-2020 counties (FALSE; default NI) or by post-2020 counties (TRUE)
counties2020 <- FALSE

# Recruitment per adult pair (FALSE; default) or per adult female (TRUE)
R_perF <- FALSE

# Drop observations of juveniles with no adults present
R_parent_drop0 <- TRUE

# Aggregation level for reproduction data
# NOTE: if this is not defined, will default to group level
sumR.Level <- "line" # Summing at the line level

# Random effects shared across areas (default TRUE)
shareRE <- FALSE

# Estimate time variation in survival (default FALSE)
survVarT <- FALSE

# Include rodent covariate on reproduction (default FALSE)
fitRodentCov <- FALSE


# DOWNLOAD/FETCH DATA #
#---------------------#

if(downloadData) {
  # (Re-)download data from GBIF
  Rype_arkiv <- downloadLN(datasets = c("Fjellstyrene", "Statskog", "FeFo"), versions = c(1.7, 1.8, 1.12), save = FALSE)
} else {
  # Use existing data sets if available
  # NOTE: Warnings about existing temporary downloads in AppData have been suppressed
  Rype_arkiv <- list()
  data.exists <- 0
  
  if(file.exists("data/FeFo_arkiv.zip")) {
    Rype_arkiv[[1]] <- suppressWarnings(LivingNorwayR::initializeDwCArchive("data/FeFo_arkiv.zip"))
    data.exists <- data.exists + 1
  }
  if(file.exists("data/Statskog_arkiv.zip")) {
    Rype_arkiv[[2]] <- suppressWarnings(LivingNorwayR::initializeDwCArchive("data/Statskog_arkiv.zip"))
    data.exists <- data.exists + 1
  }
  if(file.exists("data/Fjellstyrene_arkiv.zip")) {
    Rype_arkiv[[3]] <- suppressWarnings(LivingNorwayR::initializeDwCArchive("data/Fjellstyrene_arkiv.zip"))
    data.exists <- data.exists + 1
  }
  
  if(identical(data.exists, as.numeric(0))) {
    stop("No data sets found in /data/ folder. Set 'downloadData' to TRUE to (re-)download the data sets from GBIF.")
  } else {
    cat(paste("  ", data.exists, " data archives found and imported.", "\n", sep = ""))
  }
}


# WRANGLE LINE TRANSECT DATA #
#----------------------------#

## Set localities/areas/counties and time period of interest
areas <- listAreas()

if(counties2020) {
  counties <- listCounties2020()
} else {
  counties <- listCounties()
}

minYear <- 2007
maxYear <- 2021

## List duplicate transects to remove
duplTransects <- listDuplTransects()

## Extract transect and observational data from DwC archive
# NOTE: for county aggregation, the data must first be wrangled using area
# aggregation/filtering
LT_data_orig <- wrangleData_LineTrans(DwC_archive_list = Rype_arkiv, 
                                 duplTransects = duplTransects,
                                 #localities = localities,
                                 areas = areas,
                                 areaAggregation = TRUE,
                                 minYear = minYear,
                                 maxYear = maxYear)

## Assign transect lines and observations to counties by geographical coordinates
# NOTE: Set 'counties = NULL' to avoid filtering on listCounties()/listCounties2020()
LT_data <- assignCounty(LT_data_orig,
                        counties = counties,
                        counties2020 = counties2020)


# WRANGLE KNOWN FATE CMR DATA #
#-----------------------------#

if(survVarT) {
  ## Read in and reformat CMR data
  d_cmr <- wrangleData_CMR(minYear = minYear)
}


# WRANGLE RODENT DATA #
#---------------------#

if(fitRodentCov) {
  ## Load and reformat rodent data
  d_rodent <- wrangleData_Rodent(duplTransects = duplTransects,
                                 #localities = localities,
                                 areas = areas,
                                 areaAggregation = TRUE,
                                 minYear = minYear, maxYear = maxYear)
}

# PREPARE INPUT DATA FOR INTEGRATED MODEL #
#-----------------------------------------#

## Reformat data into vector/array list for analysis with Nimble
# NOTE: set 'd_cmr' and/or 'd_rodent' to NULL if excluding covariates
#       and time variation in survival
input_data <- prepareInputData(d_trans = LT_data$d_trans, 
                               d_obs = LT_data$d_obs,
                               d_cmr = NULL, #d_cmr
                               d_rodent = NULL, #d_rodent
                               #localities = localities, 
                               #areas = areas,
                               counties = counties,
                               areaAggregation = FALSE,
                               countyAggregation = TRUE,
                               excl_neverObs = TRUE,
                               R_perF = R_perF,
                               R_parent_drop0 = R_parent_drop0,
                               sumR.Level = sumR.Level,
                               dataVSconstants = TRUE,
                               save = TRUE)


# MODEL SETUP #
#-------------#

## Determine correct code path
code.path <- selectCodePath(shareRE = shareRE,
                            survVarT = survVarT)

## Setup for model using nimbleDistance::dHN
model_setup <- setupModel(modelCode.path = code.path,
                          customDist = TRUE,
                          nim.data = input_data$nim.data,
                          nim.constants = input_data$nim.constants,
                          R_perF = R_perF,
                          shareRE = shareRE,
                          survVarT = survVarT,
                          fitRodentCov = fitRodentCov,
                          testRun = TRUE,
                          nchains = 3,
                          initVals.seed = 0)

# MODEL (TEST) RUN #
#------------------#
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
                       setSeed = 0)
Sys.time() - t.start

saveRDS(IDSM.out, file = 'rypeIDSM_dHN_multiArea_sepRE_testRun.rds')


# TIDY UP POSTERIOR SAMPLES #
#---------------------------#
IDSM.out.tidy <- tidySamples(IDSM.out = IDSM.out, save = TRUE)


# PREPARE OUTPUT DATA FOR NI #
#----------------------------#
# Format population densities from posterior samples
prepared.output <- prepareOutputNI(
  mcmc.out = IDSM.out.tidy,
  N_areas = input_data$nim.constant$N_areas,
  area_names = input_data$nim.constant$area_names,
  N_sites = input_data$nim.constant$N_sites,
  min_years = input_data$nim.constant$min_years,
  max_years = input_data$nim.constant$max_years,
  minYear = minYear,
  maxYear = maxYear
)

# Calculate abundance (in 1000 ind.) based on estimated densities and habitat information
# NOTE: This part depends on "ptarmigan_habitat.rds", which is the output from the
#       ptarmigan habitat model by Kvasnes et al. 2018 BMC Ecology (with updated data).
#       If the rds file is absent, the below function will utilize the provided data table
#       (data/county_habitat_data.csv) with values from the last run of the model.
output.data <- estimatePtarmiganAbundance(output = prepared.output)


# WRANGLE OUTPUT DATA TO FIT NI DATABASE #
#----------------------------------------#
#devtools::install_github("NINAnor/NIcalc", build_vignettes = T)

# retrieve the currently stored values for Ptarmigan from the Nature Index (NI) database
species <- c("Lagopus lagopus")
indicators <- c("Lirype")
currentPtarmTable <- downloadData_NIdb(species, indicators, save = T, save_path = "data")

# update the values with the results from the population density model for Ptarmigan
# time interval to update can be specified using the 'min_year' and 'max_year' arguments
newPtarmTable <- updateNItable(model.est = output.data,
                               cur.table = currentPtarmTable,
                               save = T,
                               save_path = "data",
                               min_year = NULL,
                               max_year = NULL)

# upload the updated table and overwrite existing data -- NB! Not tested yet!
#uploadData_NIdb(species, data_path = "data")


# OPTIONAL: MCMC TRACE PLOTS #
#----------------------------#
# NOTE: Takes about 30+ minutes to generate PDF's

plotMCMCTraces(mcmc.out = IDSM.out.tidy,
               fitRodentCov = fitRodentCov)


# OPTIONAL: TIME SERIES PLOTS #
#-----------------------------#

# NOTE: VitalRates = TRUE will only plot reproduction if survVarT is set to FALSE
plotTimeSeries(mcmc.out = IDSM.out.tidy, 
               N_areas = input_data$nim.constant$N_areas, 
               area_names = input_data$nim.constant$area_names, 
               N_sites = input_data$nim.constant$N_sites, 
               min_years = input_data$nim.constant$min_years, 
               max_years = input_data$nim.constant$max_years, 
               minYear = minYear,
               maxYear = maxYear,
               VitalRates = TRUE,
               DetectParams = TRUE,
               Densities = TRUE,
               survVarT = survVarT)


# OPTIONAL: PLOT VITAL RATE POSTERIORS #
#--------------------------------------#

# plotPosteriorDens_VR(mcmc.out = IDSM.out.tidy,
#                      N_areas = input_data$nim.constant$N_areas, 
#                      area_names = input_data$nim.constant$area_names, 
#                      N_years = input_data$nim.constant$N_years,
#                      minYear = minYear,
#                      survAreaIdx = input_data$nim.constants$SurvAreaIdx,
#                      survVarT = survVarT,
#                      fitRodentCov = fitRodentCov) 


# OPTIONAL: PLOT COVARIATE PREDICTIONS #
#--------------------------------------#

# if(fitRodentCov){
#   plotCovPrediction(mcmc.out = IDSM.out.tidy,
#                     effectParam = "betaR.R",
#                     covName = "Rodent occupancy",
#                     minCov = 0, 
#                     maxCov = 1,
#                     N_areas = input_data$nim.constant$N_areas, 
#                     area_names = input_data$nim.constant$area_names,
#                     fitRodentCov = fitRodentCov)
# }


# OPTIONAL: MAP PLOTS #
#---------------------#

## Make map of Norwegian municipalities ("fylke")
# NorwayMunic.map <- setupMap_NorwayMunic(shp.path = "data/Kommuner_2018_WGS84/Kommuner_2018_WGS84.shp",
#                                         d_trans = LT_data$d_trans,
#                                         areas = areas, areaAggregation = TRUE)

## Plot population growth rate, density, and vital rates on map
# plotMaps(mcmc.out = IDSM.out.tidy, 
#          mapNM = NorwayMunic.map,
#          N_areas = input_data$nim.constant$N_areas, 
#          area_names = input_data$nim.constant$area_names, 
#          N_sites = input_data$nim.constant$N_sites, 
#          min_years = input_data$nim.constant$min_years, 
#          max_years = input_data$nim.constant$max_years, 
#          minYear = minYear, maxYear = maxYear,
#          fitRodentCov = fitRodentCov)


# OPTIONAL: MODEL COMPARISON (PLOTS) #
#------------------------------------#

# modelComp <- plotModelComparison(modelPaths = c("rypeIDSM_realData_Lierne.rds", 
#                                                 "rypeIDSM_dHN_realData_Lierne.rds"), 
#                                  modelChars = c("Zeroes trick", "dHN"), 
#                                  N_sites = 58, N_years = 6,
#                                  plotPath = "Plots/ModelCompTest",
#                                  returnData = FALSE)

