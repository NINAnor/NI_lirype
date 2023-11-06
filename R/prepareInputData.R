#' Prepare line transect, known fate CMR, and rodent covariate data for
#' integrated analysis
#'
#' @param d_trans tibble containing information on transects (events). Output of
#' wrangleData_LineTrans(). 
#' @param d_obs tibble containing information on observations made along
#' transects. Output of wrangleData_LineTrans().
#' @param d_cmr OPTIONAL list with 2 elements. Surv1 and Surv2 are matrices of
#' individuals released (column 1) and known to have survived (column 2) in each
#' year (row) for season 1 and 2, respectively. Output of wrangleData_CMR().
#' Set to 'NULL' to exclude time variation in survival.
#' @param d_rodent OPTIONAL matrix containing the average number of transects
#' with rodent observations per area and year. Output of wrangleData_Rodent().
#' Set to 'NULL' to exclude this covariate.
#' @param localities vector of strings listing localities to consider. Either
#' localities, areas or counties must be provided. 
#' @param areas vector of strings listing areas to consider. Either localities,
#' areas or counties must be provided.
#' @param counties vector of strings listing counties to consider. Either
#' localities, areas or counties must be provided.
#' @param areaAggregation logical. If TRUE, areas are used as smallest spatial
#' unit. If FALSE, and countyAggregation is FALSE, locations (within areas) are
#' used as smallest spatial unit.
#' @param countyAggregation logical. If TRUE, data are aggregated to regional
#' level (counties). If FALSE, and areaAggregation is FALSE, locations (within
#' areas) are used as the smallest spatial unit.
#' @param excl_neverObs logical. If TRUE (default), transects where ptarmigans
#' were never observed are excluded. If FALSE, all transects are included.
#' @param R_perF logical. If TRUE, treats recruitment rate as juvenile per adult
#' female. If FALSE, treats recruitment rate as juvenile per adult (sum of both
#' sexes).
#' @param R_parent_drop0 logical. If TRUE (default), removes observations of
#' juveniles without adults from the recruitment data. If FALSE, sets 1 as the
#' number of adults/adult females when none are observed.
#' @param sumR.Level character string. Default ("group") summarises
#' reproduction/recruitment data at the group/observation level. Setting to
#' "line" summarises data at the transect line level instead.
#' @param dataVSconstants logical. If TRUE (default) returns a list of 2 lists
#' containing data and constants for analysis with Nimble. If FALSE, returns a
#' list containing all data and constants. 
#' @param save logical. If TRUE (default) saves the prepared data in the working 
#' directory as a .rds file.
#'
#' @return A list or list of lists, depending on argument `dataVSconstants` 
#' (see above).
#' @export
#'
#' @examples

prepareInputData <- function(d_trans,
                             d_obs,
                             d_cmr,
                             d_rodent,
                             localities = NULL,
                             areas = NULL,
                             counties = NULL,
                             areaAggregation,
                             countyAggregation,
                             excl_neverObs = TRUE,
                             R_perF,
                             R_parent_drop0,
                             sumR.Level = "group",
                             dataVSconstants = TRUE,
                             save = TRUE) {
  
  # Multi-area setup #
  #------------------#
  
  ## Assignment of spatial units
  if(countyAggregation) {
    sUnits <- counties
  } else if(areaAggregation) {
    sUnits <- areas
  } else {
    sUnits <- localities
  }
  N_sUnits <- length(sUnits)
  
  ## Rename appropriate column in line transect data to reflect level of spatial aggregation
  if(countyAggregation) {
    colnames(d_trans)[which(colnames(d_trans) == "county")] <- "spatialUnit"
    colnames(d_obs)[which(colnames(d_obs) == "county")] <- "spatialUnit"
  } else if(areaAggregation){
    colnames(d_trans)[which(colnames(d_trans) == "verbatimLocality")] <- "spatialUnit"
    colnames(d_obs)[which(colnames(d_obs) == "verbatimLocality")] <- "spatialUnit"
  } else {
    colnames(d_trans)[which(colnames(d_trans) == "locality")] <- "spatialUnit"
    colnames(d_obs)[which(colnames(d_obs) == "locality")] <- "spatialUnit"
  }
  
  ## If desired: remove all transects on which willow ptarmigans were never encountered --> of all time, or per year?
  if(excl_neverObs){
    
    #This removes transects without observations in a given year, not only
    #transects where Willow ptarmigan was never encountered
    #d_trans <- dplyr::semi_join(d_trans, d_obs, by = c("locationID", "Year"))
    
    #By removing the 'Year' parameter, only transects without any observations
    #of Willow ptarmigan across all years are removed
    d_trans <- dplyr::semi_join(d_trans, d_obs, by = c("locationID"))
  }

  ## Variables shared across areas
  # Number of age classes
  N_ageC <- 2
  
  # (Number of) years
  range_yearsTot <- min(d_trans$Year):max(d_trans$Year)
  idx_yearsTot <- range_yearsTot - min(d_trans$Year) + 1
  N_yearsTot <- length(range_yearsTot)
  
  # Truncation distance
  W <- 200
  
  ## Count sites per area (for defining correct array dimensions)
  site_count <- d_trans %>%                             
    dplyr::group_by(spatialUnit) %>%
    dplyr::summarise(count = n_distinct(locationID))
  
  if(length(which(site_count$count == 0)) > 0){
    stop("One or more specified areas contain no willow ptarmigan observations")
  }
  
  ## Count observations per spatial unit
  obs_count <- d_obs %>% 
    dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
    dplyr::group_by(spatialUnit) %>%
    dplyr::summarise(count = n())
    
  ## Set up arrays for spatial unit specific data
  N_sites <- N_years <- min_years <- max_years <- N_obs <- N_sumR_obs <- rep(NA, N_sUnits)
  
  A <- matrix(NA, nrow = N_sUnits, ncol = N_yearsTot)
  y <- Year_obs <- zeros_dist <- matrix(NA, nrow = N_sUnits, ncol = max(obs_count$count))
  sumR_obs <- sumR_obs_year <- sumAd_obs<- matrix(NA, nrow = N_sUnits, ncol = max(obs_count$count))
  
  L <- N_line_year <- N_J_line_year <- N_A_line_year <- array(0, dim = c(N_sUnits, max(site_count$count), N_yearsTot))
  
  N_a_line_year <- array(0, dim = c(N_sUnits, N_ageC, max(site_count$count), N_yearsTot))
  
  
  for(x in 1:N_sUnits){
    
    ## Subset data (specific area)
    # if(areaAggregation){
    #   d_trans_sub <- subset(d_trans, spatialUnit == areas[x])
    #   d_obs_sub <- subset(d_obs, spatialUnit == areas[x])
    # }else{
    #   d_trans_sub <- subset(d_trans, spatialUnit == localities[x])
    #   d_obs_sub <- subset(d_obs, spatialUnit == localities[x])
    # }
    if(!(sUnits[x] %in% d_trans$spatialUnit)){
      stop(paste0("Spatial unit ", sUnits[x], " (index ", x, ") is not in the data."))
    }
    
    d_trans_sub <- subset(d_trans, spatialUnit == sUnits[x])
    d_obs_sub <- subset(d_obs, spatialUnit == sUnits[x])
    
    # Constants #
    #-----------#
    
    ## Numbers of sites
    N_sites[x] <- n_distinct(d_trans_sub$locationID)
    
    ## (Number of) years with data
    years_mon <- sort(unique(d_trans_sub$Year)) - min(range_yearsTot) + 1 
    min_years[x] <- min(years_mon)
    max_years[x] <- max(years_mon)

    # Transect characteristics #
    #--------------------------#
    
    ## Tansect lengths
    TransLen <- d_trans_sub %>% 
      dplyr::select(locationID, Year, sampleSizeValue) %>%
      reshape2::dcast(locationID~Year, value.var = "sampleSizeValue", sum) %>%
      arrange(locationID)
    
    TransLen_mat <- TransLen %>% select(-locationID) %>% as.matrix()
    colnames(TransLen_mat) <- as.numeric(colnames(TransLen_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        L[x, 1:N_sites[x], t] <- TransLen_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
 
    
    ## Total covered area 
    A[x, 1:N_yearsTot] <- colSums(L[x,,])*W*2
    
    
    # Observation distance from transect #
    #------------------------------------#
 
    
    temp_dist <- d_obs_sub %>% 
      dplyr::filter(between(DistanceToTransectLine, -0.1, W)) %>% 
      dplyr::select(Year, DistanceToTransectLine) %>%
      dplyr::mutate(Year2 = Year - (minYear) + 1)
    
    ## Number of observations
    N_obs[x] <- length(temp_dist$DistanceToTransectLine)
    
    ## Distance to transect line
    y[x, 1:N_obs[x]] <- temp_dist$DistanceToTransectLine
    
    ## Observation year
    Year_obs[x, 1:N_obs[x]] <- temp_dist$Year2
    
    ## Vector of 0's of same length as y
    zeros_dist[x, 1:N_obs[x]] <- rep(0, N_obs[x])
    
    
    # Number of birds/line (pooled age classes) #
    #-------------------------------------------#
    
    temp <- TransLen %>% select(locationID)

    TaksObs <- d_obs_sub %>% 
      filter(between(DistanceToTransectLine, -0.1, W)) %>% # 
      dplyr::mutate(cs = unknownJuvenile+unknownunknown+FemaleAdult+MaleAdult) %>%
      reshape2::dcast(locationID~Year, value.var = "cs", sum) %>%
      dplyr::right_join(., temp, by = c("locationID" = "locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
  
    TaksObs_mat <- TaksObs %>% dplyr::select(-locationID) %>% as.matrix() 
    colnames(TaksObs_mat) <- as.numeric(colnames(TaksObs_mat)) - min(range_yearsTot) + 1
    
    # filtering observations on DistanceToTransectLine will in some cases remove
    # all observations for a transect line within 1 year. This recalculates the
    # number of years with samples
    years_mon <- as.numeric(colnames(TaksObs_mat))
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_line_year[x, 1:N_sites[x], t] <- TaksObs_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
    
    # Number of birds/line (by age class) #
    #-------------------------------------#
    
    ## Juveniles (& unknowns)
    TaksObs_J <- d_obs_sub %>% 
      filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = unknownJuvenile + unknownunknown) %>%
      reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
      dplyr::right_join(., temp, by=c("locationID"="locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
    
    TaksObs_J_mat <- TaksObs_J %>% dplyr::select(-locationID) %>% as.matrix() 
    colnames(TaksObs_J_mat) <- as.numeric(colnames(TaksObs_J_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_J_line_year[x, 1:N_sites[x], t] <- TaksObs_J_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
    
    ## Adults 
    TaksObs_A <- d_obs_sub %>% 
      filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(cs = FemaleAdult + MaleAdult) %>%
      reshape2::dcast(locationID~Year, value.var="cs", sum) %>%
      dplyr::right_join(., temp, by=c("locationID"="locationID")) %>%
      replace(., is.na(.), 0) %>%
      dplyr::arrange(locationID)
    
    TaksObs_A_mat <- TaksObs_A %>% dplyr::select(-locationID) %>% as.matrix() 
    colnames(TaksObs_A_mat) <- as.numeric(colnames(TaksObs_A_mat)) - min(range_yearsTot) + 1
    
    for(t in 1:N_yearsTot){
      
      if(t %in% years_mon){
        N_A_line_year[x, 1:N_sites[x], t] <- TaksObs_A_mat[1:N_sites[x], which(years_mon == t)]
      }
    }
    
    ## Check and combine in array
    if(!(all(N_J_line_year + N_A_line_year == N_line_year))){
      warning("Number of observed adults and juveniles does not add up correctly. Double-check data.")
    }
    
    N_a_line_year[x, 1,,] <- N_J_line_year[x,,]
    N_a_line_year[x, 2,,] <- N_A_line_year[x,,]
    
    
    # Recruitment #
    #-------------#
    
    # NOTE: For consistency with the population model, we need to define R as the number of recruits
    # per adult (female)
    
    # Reformat data
    temp_Rec <- d_obs_sub %>% filter(between(DistanceToTransectLine, -0.1, W)) %>%
      dplyr::mutate(sumR = unknownJuvenile + unknownunknown,
                    sumAd = MaleAdult + FemaleAdult,
                    sumAdF = FemaleAdult) %>%
      dplyr::mutate(Year2 = Year - (min(Year)) + 1)
  
    # Optional: summarise data at line level (per year)
    if(sumR.Level == "line"){
      temp_Rec <- temp_Rec %>%
        dplyr::group_by(locationID, Year, Year2) %>%
        dplyr::summarise(sumR = sum(sumR),
                         sumAd = sum(sumAd),
                         sumAdF = sum(sumAdF), .groups = "keep")
    }
    
    # Set which adult measure to use
    if(R_perF){
      temp_Rec$sumAd_use <- temp_Rec$sumAdF
    }else{
      temp_Rec$sumAd_use <- temp_Rec$sumAd
    }
    
    # Deal with instances of 0 adults observed
    if(R_parent_drop0){ # --> Drop all cases of 0 adults observed
      temp_Rec <- temp_Rec %>%
        dplyr::filter(sumAd_use != 0)
      
    }else{ # --> Drop only cases when neither adult (females) nor juveniles were observed and add +1 to adults otherwise
      
      temp_Rec <- temp_Rec %>%
        dplyr::filter(sumAd_use + sumR != 0) %>%
        dplyr::mutate(sumAd_use = ifelse(sumAd_use == 0, 1, sumAd_use))
    }

    # Extract relevant data vectors
    N_sumR_obs[x] <- length(temp_Rec$sumR)
    sumR_obs[x, 1:N_sumR_obs[x]] <- temp_Rec$sumR
    sumAd_obs[x, 1:N_sumR_obs[x]] <- temp_Rec$sumAd_use
    sumR_obs_year[x, 1:N_sumR_obs[x]] <- temp_Rec$Year2
  }
  
  ## Drop excess columns in recruitment data
  sumR_obs <- sumR_obs[1:N_sUnits, 1:max(N_sumR_obs), drop = FALSE]
  sumAd_obs <- sumAd_obs[1:N_sUnits, 1:max(N_sumR_obs), drop = FALSE]
  sumR_obs_year <- sumR_obs_year[1:N_sUnits, 1:max(N_sumR_obs), drop = FALSE]
  
  
  # Data assembly #
  #---------------#
  
  ## If including survival estimates, set spatial index for CMR data
  if(!is.null(d_cmr)) {
    if(countyAggregation) {
      SurvAreaIdx <- which(sUnits == d_cmr$county_names)
    } else if(areaAggregation) {
      SurvAreaIdx <- which(sUnits == d_cmr$area_names)
    } else {
      SurvAreaIdx <- which(sUnits == d_cmr$locality_names)
    }
    
    if(length(SurvAreaIdx) == 0){
      stop("No overlap in areas for line transect and survival data.\n
           The present implementation of the model requires including\n
           line transect data from Lierne.")
    }
  }
  
  ## Add dummy dimensions if running for only one spatial unit
  if(N_sUnits == 1){
    N_sites <- c(N_sites, NA)
  }
  
  ## Assembling all data in a list
  input.data <- list(
    sumR_obs = sumR_obs, # Observed numbers of recruits
    sumAd_obs = sumAd_obs, # Observed numbers of adults/adult females
    sumR_obs_year = sumR_obs_year, # Year of observed numbers of recruits
    N_sumR_obs = N_sumR_obs, # Total number of observations of numbers of recruits
    
    y = y, # Distance to transect line for each individual observation
    zeros_dist = zeros_dist, # Vector of 0's of same length as y
    Year_obs = Year_obs, # Year of each observation
    N_obs = N_obs, # Total number of observations
    
    N_line_year = N_line_year, # Number of birds observed per site per year
    N_a_line_year = N_a_line_year, # Number of birds observed per ageclass per site per year
    L = L, # Transect length per site and year
    
    N_years = N_yearsTot, # Max number of years with data
    min_years = min_years, # Earliest year (index) of monitoring per area
    max_years = max_years, # Last year (index) of monitoring per area
    N_sites = N_sites, # Total number of monitored sites per area
    
    A = A, # Total covered area per year
    W = W, # Truncation distance
    N_ageC = N_ageC, # Number of age classes
    
    N_areas = N_sUnits, # Number of areas (integer)
    area_names = sUnits # Area names (character vector)
  )
  
  if (!is.null(d_cmr)) {
    input.data$Survs1 <- d_cmr$Survs1 # Season 1 releases & survivors (area 1)
    input.data$Survs2 <- d_cmr$Survs2 # Season 2 releases & survivors (area 1)
    input.data$SurvAreaIdx <- SurvAreaIdx
    input.data$year_Survs <- d_cmr$year_Survs # Years (indices) of telemetry data
    input.data$N_years_RT <- length(d_cmr$year_Survs)
  }
  
  if (!is.null(d_rodent)) {
    input.data$RodentOcc <- d_rodent # Rodent covariate (matrix)
  }
  
  ## Assembling Nimble data
  nim.data <- list(
    sumR_obs = input.data$sumR_obs,
    sumAd_obs = sumAd_obs,
    y = input.data$y,
    zeros.dist = input.data$zeros_dist,
    L = input.data$L,
    N_line_year = input.data$N_line_year,
    N_a_line_year = input.data$N_a_line_year,
    A = input.data$A
  )
  
  if (!is.null(d_cmr)) {
    nim.data$Survs1 <- input.data$Survs1
    nim.data$Survs2 <- input.data$Survs2
  }
  
  if (!is.null(d_rodent)) {
    nim.data$RodentOcc <- input.data$RodentOcc
  }
  
  ## Assembling Nimble constants
  nim.constants <- list(
    N_years = input.data$N_years,
    min_years = input.data$min_years,
    max_years = input.data$max_years,
    W = input.data$W,
    N_obs = input.data$N_obs,
    Year_obs = input.data$Year_obs,
    N_sites = input.data$N_sites,
    N_ageC = N_ageC,
    N_areas = input.data$N_areas,
    area_names = input.data$area_names,
    sumR_obs_year = input.data$sumR_obs_year,
    N_sumR_obs = input.data$N_sumR_obs
  )
  
  if (!is.null(d_cmr)) {
    nim.constants$SurvAreaIdx <- input.data$SurvAreaIdx
    nim.constants$year_Survs <- input.data$year_Survs
    nim.constants$N_years_RT <- input.data$N_years_RT
  }

  ## Make final data list to return
  if (dataVSconstants) {
    rype.data <- list(nim.data = nim.data,
                      nim.constants = nim.constants)
  } else {
    rype.data <- input.data
  }
  
  ## Optional: save data as .rds
  if (save) {
    saveRDS(rype.data, file = "RypeData_forIM.rds")
  }
  
  ## Return data
  return(rype.data)
}

