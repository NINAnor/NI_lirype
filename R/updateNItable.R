#' Prepare updated indicator data for the Nature Index database
#'
#' This function wrangles the currently stored values for Ptarmigan from the
#' Nature Index (NI) database, and updates the values according to the output
#' from the population model. Note that the time period for updates can be
#' specified using the `min_year` and `max_year` arguments. Values outside this
#' range will be copied from the existing indicator values. If the time period
#' is not specified, all indicator values from and including 2007 will be updated.
#' 
#' @param output.data a data frame with estimated abundances, lower 25% quantile
#' and upper 75% quantile for ptarmigan per county and year. Values are in 1000
#' individuals.
#' @param currentPtarmTable a list with the downloaded ptarmigan table from
#' the Nature Index (NI) database.
#' @param save logical, default = `FALSE`. If `TRUE`, the updated indicator 
#' data is saved in the workspace in a file named `updatedIndicatorData.rds`. 
#' @param save_path character string indicating the directory into which to save
#' the updated indicator data.
#' @param min_year numeric, sets the starting year for updates to the indicator
#' values. If `max_year` is not set, all values from and including `min_year`
#' will be updated.
#' @param max_year numeric, sets the end year for updates to the indicator
#' values. If `min_year` is not set, all the values up to and including
#' `max_year` will be updated.
#'
#' @return A list containing 1) a data frame of the values of the indicator
#' and 2) a list of distribution objects quantifying the uncertainty of each
#' indicator value estimate. Only 1) is updated here, whereas 2) is attached
#' from the original download for compatibility.
#' @export
#'
#' @examples
#' 
updateNItable <- function(
    model.est = output.data,
    cur.table = currentPtarmTable,
    save = FALSE,
    save_path,
    min_year = NULL,
    max_year = NULL) {
  
  ind.values <- cur.table$Lirype$indicatorValues

  # filter model estimates by min_year, if specified
  if(!is.null(min_year)) {
    if(is.numeric(min_year) & dplyr::between(min_year, min(model.est$Year), max(model.est$Year))) {
      model.est <- model.est %>% dplyr::filter(Year >= min_year)
    } else {
      stop(paste0("`min_year` value must be between ", min(model.est$Year), " and ", max(model.est$Year)))
    }
  }
  
  # filter model estimates by max_year, if specified
  if(!is.null(max_year)) {
    if(is.numeric(max_year) & dplyr::between(max_year, min(model.est$Year), max(model.est$Year))) {
      model.est <- model.est %>% dplyr::filter(Year <= max_year)
    } else {
      stop(paste0("`max_year` value must be between ", min(model.est$Year), " and ", max(model.est$Year)))
    }
  }
  
  # conform model estimates table with indicator values table
  colnames(model.est) <- c("areaName", "yearName", "verdi", "nedre_Kvartil", "ovre_Kvartil")
  model.est$yearName <- as.character(model.est$yearName)  
  
  # create reference tables for area ID's/names and year ID's/names
  area.table <- dplyr::distinct(ind.values[,c("areaId","areaName")])
  year.table <- data.frame(sortID = seq(1,length(unique(ind.values$yearName)),1),
                           yearName = sort(unique(ind.values$yearName), decreasing = T))
  year.table <- dplyr::left_join(year.table, ind.values[,c("yearName","yearId")], by = "yearName")
  year.table <- year.table[!is.na(year.table$yearId),]
  year.table <- dplyr::distinct(year.table)
  
  # add year sortID to keep current table order
  ind.values <- dplyr::left_join(ind.values, year.table, by = "yearName")
  
  # replace the rows in the current ptarmigan table matching area and year in
  # the model output table
  ind.values <- ind.values %>% 
    dplyr::anti_join(model.est, by = c("areaName", "yearName")) %>%
    dplyr::bind_rows(model.est)
  
  # match areaID/name and yearID/name to newly added rows from ref.tables
  ind.values <- dplyr::left_join(ind.values, area.table, by = "areaName")
  ind.values <- dplyr::left_join(ind.values, year.table, by = "yearName")

  # sort table by area and year
  ind.values <- ind.values[order(ind.values$areaId.y, ind.values$sortID.y, decreasing = F), ]
  
  # add missing values for indicatorId, indicatorName and unitOfMeasurement
  ind.values$indicatorId[is.na(ind.values$indicatorId)] <- 109
  ind.values$indicatorName[is.na(ind.values$indicatorName)] <- "Lirype"
  ind.values$unitOfMeasurement[is.na(ind.values$unitOfMeasurement)] <- "thousand ind"
  
  # yet to deal with: datatypeId and datatypeName !!
  
  # sort and rename columns
  ind.values <- ind.values %>% dplyr::select("indicatorId", "indicatorName", "areaId.y",
                                             "areaName", "yearId", "yearName",
                                             "verdi", "nedre_Kvartil", "ovre_Kvartil",
                                             "datatypeId", "datatypeName", "unitOfMeasurement",
                                             "customDistributionUUID", "distributionName",
                                             "distributionId", "distParam1", "distParam2")
  
  # fix column names
  colnames(ind.values)[3] <- "areaId"
  
  # add updated indicator values to the original list structure
  cur.table$Lirype$indicatorValues <- ind.values
  
  # save the updated indicator data
  if(save){
    saveRDS(cur.table, file = paste0(save_path, "/updatedIndicatorData.rds"))
  }
  
  # return
  return(cur.table)
  
}