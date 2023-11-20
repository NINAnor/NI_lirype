#' Read in and reformat known fate capture-mark-recapture (CMR) data
#'
#' @param minYear integer. Earliest year of line transect data used.
#'
#' @return list with 2 elements. Surv1 and Surv2 are matrices of individuals 
#' released (column 1) and known to have survived (column 2) in each year (row)
#' for season 1 and season 2, respectively.
#' @export
#'
#' @examples
#'     wrangleData_CMR(minYear = minYear)

wrangleData_CMR <- function(minYear){
  
  ## Load CMR data
  CMR_data <- tibble::as_tibble(read.csv("data/CMR_Data.csv", header = T, sep = ",")) 
  
  ## Divide data between time periods
  
  # Period 1
  Survs1 <- CMR_data %>% 
    dplyr::filter(TimePeriod == 1) %>% 
    dplyr::select(-YearPeriod, -TimePeriod) %>%
    as.matrix()
  
  # Period 2
  Survs2 <- CMR_data %>% 
    dplyr::filter(TimePeriod == 2) %>% 
    dplyr::select(-YearPeriod, -TimePeriod) %>%
    as.matrix()
  
  # Time indices
  year_Survs <- unique(CMR_data$YearPeriod) - minYear + 1
  
  ## Arrange in list
  d_cmr <- list(Survs1 = Survs1, 
                Survs2 = Survs2, 
                year_Survs = year_Survs,
                locality_names = "Lierne Fjellst. Vest",
                area_names = "Lierne Fjellstyre")
  
  ## Return list
  return(d_cmr)
}