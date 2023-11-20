#' Extract and filter transect and observation data from DwC archive
#'
#' @param DwC_archive_list list of DwCArchives. Darwin Core archives containing
#' line transect data for willow ptarmigan (Lagopus lagopus) in Norway.
#' @param duplTransects vector of strings. (parent) EventIDs of transects that
#' are duplicates and need removing from the data
#' @param localities string or vector of strings. Names of localities to extract
#' data for. Either localities or areas must be provided. 
#' @param areas string or vector of strings. Names of areas to extract
#' data for. Either localities or areas must be provided.
#' @param areaAggregation logical. If TRUE, areas are used as smallest spatial
#' unit. If FALSE, locations (within areas) are used as smallest spatial unit.
#' @param minYear integer. Earliest year of data to extract.
#' @param maxYear integer. Latest year of data to extract.  
#'
#' @return list of 2 tibbles. `d_trans` contains information on transects 
#' (events). `d_obs`contains information on observations made along transects 
#' (distance to transect line, numbers of birds in each age/sex class observed,
#' etc.)
#' 
#' @export
#'
#' @examples
#'    wrangleData_LineTrans(DwC_archive_list = Rype_arkiv,
#'                          duplTransects = duplTransects,
#'                          areas = areas,
#'                          areaAggregation = TRUE,
#'                          minYear = minYear,
#'                          maxYear = maxYear)

wrangleData_LineTrans <- function(DwC_archive_list,
                                  duplTransects,
                                  localities = NULL,
                                  areas = NULL,
                                  areaAggregation,
                                  minYear,
                                  maxYear) {
    
  ## Extract relevant parts from DwC_archive
  
  # Set up lists for storing data
  Eve_all <- Occ <- Meas <- list()
  
  for(i in 1:length(DwC_archive_list)) {
    # Core table
    Core <- DwC_archive_list[[i]]$getCoreTable()
    
    # Complete event table
    # NOTE: added 'suppressWarnings' on warning message 'NAs introduced by coercion' for
    #  'as.numeric(sampleSizeValue)' (regarding missing values for transect line length)
    Eve_all[[i]] <- suppressWarnings(tibble::as_tibble(Core$exportAsDataFrame()) %>%
      dplyr::mutate(sampleSizeValue = as.numeric(sampleSizeValue)), )
    
    # Complete occurrence table
    Occ[[i]] <- tibble::as_tibble(DwC_archive_list[[i]]$getExtensionTables("occurrence")[[1]]$exportAsDataFrame())
    
    # Complete measurements and facts table
    Meas[[i]] <- tibble::as_tibble(DwC_archive_list[[i]]$getExtensionTables()[[1]]$exportAsDataFrame())
  }

  # Bind datasets together
  Eve_all <- dplyr::bind_rows(Eve_all, .id = "column_label")
  Occ <- dplyr::bind_rows(Occ, .id = "column_label")
  Meas <- dplyr::bind_rows(Meas, .id = "column_label")
  
  ## Filter event data by area and year
  if(areaAggregation){
    Eve <- Eve_all %>% 
      dplyr::mutate(eventDate = as.Date(eventDate)) %>%
      dplyr::mutate(Year = lubridate::year(eventDate)) %>%
      dplyr::filter(verbatimLocality %in% areas) %>%
      dplyr::filter(dplyr::between(Year, minYear, maxYear))
  }else{
    Eve <- Eve_all %>% 
      dplyr::mutate(eventDate = as.Date(eventDate)) %>%
      dplyr::mutate(Year = lubridate::year(eventDate)) %>%
      dplyr::filter(locality %in% localities) %>%
      dplyr::filter(dplyr::between(Year, minYear, maxYear))
  }
  
  ## Assemble transect level info
  d_trans <- Eve %>% 
    dplyr::select(locationID, eventDate, eventID, modified, 
                  samplingProtocol, eventRemarks, sampleSizeValue, 
                  stateProvince, municipality, locality, 
                  verbatimLocality, locationRemarks, Year, footprintWKT) %>% #LRE: added footprintWKT for spatial filtering
    dplyr::filter(eventRemarks == "Line transect") %>%
    dplyr::mutate(locationRemarks = gsub("In the original data this is known as lineID ", '', locationRemarks)) %>%
    tidyr::separate(., col = locationRemarks, sep = ",", into = c("LineName", "locationRemarks")) %>%
    dplyr::select(-locationRemarks)
  
  ## Identify and remove transects with suspiciously short length (< 200 m) and duplicate transects
  bad_transects <- c(d_trans$eventID[which(d_trans$sampleSizeValue < 200)], duplTransects)
  d_trans <- d_trans %>%
    dplyr::filter(!(eventID %in% bad_transects))
  
  ## Double-check no duplicate transects remain
  transect_duplicates <- d_trans %>%
    dplyr::mutate(year = lubridate::year(eventDate)) %>%
    dplyr::group_by(locationID, locality, verbatimLocality, year) %>%
    dplyr::summarise(transectCount = dplyr::n(), .groups = 'keep') %>%
    dplyr::filter(transectCount > 1)
  
  if(nrow(transect_duplicates) > 0){
    stop("There are duplicate transects (> 1 transect in same location per year). Identify the duplicate transects and update R/listDuplTransects.R accordingly.")
  }
  
  ## Assemble observation level info
  # Observations: distance to transect lines
  d_obsTemp <- Eve %>% 
    dplyr::select(locationID, locality, verbatimLocality, parentEventID, eventID,
                  eventRemarks, dynamicProperties, eventDate) %>%
    dplyr::filter(eventRemarks == "Human observation" & !is.na(dynamicProperties)) %>%
    dplyr::mutate(Year = lubridate::year(eventDate)) %>%
    dplyr::left_join(., Meas, by = c("eventID" = "id")) %>%
    dplyr::select(locationID, locality, verbatimLocality, parentEventID, eventID, measurementValue, Year)
  colnames(d_obsTemp)[6] <- "DistanceToTransectLine"
  
  # Observations: remaining information (willow ptarmigan only)
  d_obs <- Occ %>% 
    dplyr::filter(!is.na(eventID)) %>%
    dplyr::select(eventID, scientificName, individualCount, sex, lifeStage) %>%
    dplyr::mutate(lifeStage = if_else(is.na(lifeStage), "unknown", lifeStage),
                  sex = if_else(is.na(sex), "unknown", sex)) %>%
    dplyr::mutate(SexStage = str_c(sex, lifeStage, sep = "")) %>%
    dplyr::select(-sex, -lifeStage) %>%
    tidyr::spread(key = "SexStage", value = "individualCount", fill = 0) %>%
    dplyr::right_join(., d_obsTemp, by = "eventID") %>%
    dplyr::filter(scientificName == "Lagopus lagopus")
  
  ## Remove observations from transects with suspiciously short length (< 200) and duplicate transects
  d_obs <- d_obs %>%
    dplyr::filter(!(parentEventID %in% bad_transects))
  
  ## Remove observations with negative or suspiciously long distances from the transect line
  d_obs <- d_obs %>%
    dplyr::filter(!DistanceToTransectLine < 0 & !DistanceToTransectLine > 1000)
  
  ## Collate and return data
  LT_data <- list(d_trans = d_trans, d_obs = d_obs)
}
