#' List counties to be included in the analysis - pre-2020 county borders
#' (corresponding to NI database)
#' 
#' @return a character vector with UTF-8 encoding. 
#' @export
#'
#' @examples
#'     listCounties()

listCounties <- function(){
  
  counties <- c(
    #"Østfold", #No ptarmigan
    #"Akershus", #No ptarmigan
    #"Oslo", #No ptarmigan
    "Hedmark",
    "Oppland",
    "Buskerud", #Only 1 area GBIF
    #"Vestfold", #No data GBIF
    #"Telemark", #No data GBIF
    #"Aust-Agder", #No data GBIF
    "Vest-Agder", #Only 1 area GBIF
    "Rogaland", #Only 1 area GBIF
    "Hordaland", #Only 1 area GBIF
    #"Sogn og Fjordane", #No data GBIF
    #"Møre og Romsdal", #No data GBIF
    "Sør-Trøndelag",
    "Nord-Trøndelag",
    "Nordland",
    "Troms",
    "Finnmark"
  )
  
  Encoding(counties) <- "UTF-8"
  return(counties)
}
