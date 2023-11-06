## Load packages required to define the pipeline
library(targets)
library(nimble)
library(nimbleDistance)
library(tidyverse)
library(ggplot2)
library(LivingNorwayR)
library(rgdal)
library(rgeos)
library(sf)
library(reshape2)
library(maptools)
library(tmap)

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
counties2020 <- FALSE #pre- (FALSE; default NI) or post- (TRUE) 2020 counties
R_perF <- FALSE # Recruitment per adult or per adult female
R_parent_drop0 <- TRUE # Drop observations of juveniles with no adults present
sumR.Level <- "line" # Aggregation level for reproduction data (line vs. group)
shareRE <- TRUE # Random effects shared across areas
survVarT <- FALSE # Time variation in survival
fitRodentCov <- FALSE # Rodent covariate on reproduction
species <- c("Lagopus lagopus") # For downloading / uploading indicator values
indicators <- c("Lirype") # For download / uploading indicator values


## Set target-specific options such as packages.
tar_option_set(packages = c("LivingNorwayR", "tidyverse", "qs"),
               format = "qs",
               memory = "transient", 
               garbage_collection = TRUE)

## Define Targets List
list(
  tar_target(
    Rype_arkiv,
    downloadLN(datasets = c("Fjellstyrene", "Statskog", "FeFo"), versions = c(1.7, 1.8, 1.12), save = TRUE)
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
    counties,
    listCounties()
  ),
  
  tar_target(
    LT_data_orig,
    wrangleData_LineTrans(DwC_archive_list = Rype_arkiv,
                          duplTransects = duplTransects,
                          #localities = localities,
                          areas = areas,
                          areaAggregation = areaAggregation,
                          minYear = minYear,
                          maxYear = maxYear)
  ),
  
  tar_target(
    LT_data,
    assignCounty(LT_data_orig,
                 counties = counties,
                 counties2020 = counties2020)
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
                       areaAggregation = areaAggregation,
                       minYear = minYear, maxYear = maxYear)
  ),
  
  tar_target(
    input_data,
    prepareInputData(d_trans = LT_data$d_trans, 
                     d_obs = LT_data$d_obs,
                     d_cmr = NULL,
                     d_rodent = NULL,
                     #localities = localities,
                     #areas = areas,
                     areaAggregation = areaAggregation,
                     countyAggregation = TRUE,
                     excl_neverObs = TRUE,
                     R_perF = R_perF,
                     R_parent_drop0 = R_parent_drop0,
                     sumR.Level = sumR.Level,
                     dataVSconstants = TRUE,
                     save = TRUE)
  ),
  
  tar_target(
    code.path,
    selectCodePath(shareRE = shareRE,
                   survVarT = survVarT)
  ),
  
  tar_target(
    model_setup,
    setupModel(modelCode.path = "NIMBLE code/rypeIDSM_multiArea_dHN.R",
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
  ),
  
  tar_target(
    IDSM.out,
    nimbleMCMC(code = model_setup$modelCode,
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
  ),
  
  tar_target(
    IDSM.out.tidy,
    tidySamples(IDSM.out = IDSM.out,
                save = TRUE)
  ),
  
  tar_target(
    prepared.output,
    prepareOutputNI(mcmc.out = IDSM.out.tidy,
                    N_areas = input_data$nim.constant$N_areas,
                    area_names = input_data$nim.constant$area_names,
                    N_sites = input_data$nim.constant$N_sites,
                    min_years = input_data$nim.constant$min_years,
                    max_years = input_data$nim.constant$max_years,
                    minYear = minYear,
                    maxYear = maxYear)
  ),
  
  tar_target(
    output.data,
    estimatePtarmiganAbundance(output = prepared.output)
  ),
  
  tar_target(
    currentPtarmTable,
    downloadData_NIdb(species,
                      indicators,
                      save = T,
                      save_path = "data")
  ),
  
  tar_target(
    newPtarmTable,
    updateNItable(model.est = output.data,
                  cur.table = currentPtarmTable,
                  save = T,
                  save_path = "data",
                  min_year = NULL,
                  max_year = NULL)
  ),
  
  # tar_target(
  #   uploadData_NI,
  #   uploadData_NIdb(species, data_path = "data")
  # ),
  
  tar_target(
    mcmc.tracePlots,
    plotMCMCTraces(mcmc.out = IDSM.out.tidy,
                   fitRodentCov = fitRodentCov),
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
                   minYear = minYear,
                   maxYear = maxYear,
                   VitalRates = TRUE,
                   DetectParams = TRUE,
                   Densities = TRUE,
                   survVarT = survVarT),
    format = "file"
  ),
  
  # tar_target(
  #   NorwayMunic.map,
  #   setupMap_NorwayMunic(shp.path = "data/Kommuner_2018_WGS84/Kommuner_2018_WGS84.shp",
  #                       d_trans = LT_data$d_trans,
  #                       areas = areas, areaAggregation = areaAggregation)
  # ),
  
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
                      N_areas = input_data$nim.constant$N_areas, 
                      area_names = input_data$nim.constant$area_names,
                      fitRodentCov = fitRodentCov)
  ),
  
  tar_target(
    mapPlots,
    plotMaps(mcmc.out = IDSM.out.tidy, 
             mapNM = NorwayMunic.map,
             N_areas = input_data$nim.constant$N_areas, 
             area_names = input_data$nim.constant$area_names, 
             N_sites = input_data$nim.constant$N_sites, 
             min_years = input_data$nim.constant$min_years, 
             max_years = input_data$nim.constant$max_years, 
             minYear = minYear, maxYear = maxYear,
             fitRodentCov = fitRodentCov)
  )
)

# Test by running tar_manifest(fields = all_of("command")) and tar_visnetwork() in the console

# Run workflow using tar_make() in the console
# Check with tar_network() in the console

# We can then use tar_read() and tar_load() to inspect and work with results
