#' Calculate the number of ptarmigan (in 1000) per spatial unit, based on
#' model estimates (number of ptarmigan per km2) and available ptarmigan habitat
#' classified as medium suitability or above (classes 4-6).
#'
#' @param output a list of model outputs (see 'prepareOutputNI.R')
#' @param rds path to ptarmigan_habitat.rds, if available
#' @param csv path to county_habitat_data.csv, if available
#'
#' @return a data frame with the estimated number of individuals (in 1000's) per
#' spatial unit per year
#' @export
#'
#' @examples

estimatePtarmiganAbundance <- function(output = prepared.output,
                                       rds = "data/ptarmigan_habitat.rds",
                                       csv = "data/county_habitat_data.csv") {
  
  if(file.exists(rds)) {
    #ptarmHabitat
    load(rds)
  } else if (file.exists(csv)) {
    ptarmHabitat <- read.csv2(file = csv, header = T, sep = ";", dec = ",", stringsAsFactors = F)
  } else {
    stop("Ptarmigan habitat data not found.")
  }
  
  # join by county
  annualArea <- output$annual.area
  ptarmHabitat <- ptarmHabitat %>% dplyr::select("countyName", "sqkm_4.6")
  annualArea <- dplyr::left_join(annualArea, ptarmHabitat, by = c("Area" = "countyName"))
  
  # calculate number of individuals in total (in 1000's)
  annualArea$estimate <- round((annualArea$Median * annualArea$sqkm_4.6) / 1000, 3)
  annualArea$lower25 <- round((annualArea$lCI * annualArea$sqkm_4.6) / 1000, 3)
  annualArea$upper75 <- round((annualArea$uCI * annualArea$sqkm_4.6) / 1000, 3)
  
  # clean and return
  annualArea <- annualArea %>% dplyr::select("Area", "Year", "estimate", "lower25", "upper75")

  # OPTIONAL: plot time series of no. of individuals per spatial unit  
  # ggplot(data = annualArea, aes(x = Year, y = estimate, group = Area)) +
  #   geom_line(aes(col = Area)) +
  #   labs(x="year", y = "No. individuals (1000)") +
  #   scale_x_continuous(breaks=seq(2007, 2021, 1))

  return(annualArea)  
}