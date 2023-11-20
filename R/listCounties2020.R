#' List counties to be included in the analysis - post-2020 county borders
#'
#' @return a character vector with UTF-8 encoding. 
#' @export
#'
#' @examples
#'     listCounties2020()

listCounties2020 <- function(){
  
  counties <- c(
    #"Agder",
    #"Rogaland",
    #"Vestland",
    #"Viken",
    #"Møre og Romsdal",
    #"Vestfold og Telemark",
    "Innlandet",
    "Nordland",
    "Troms og Finnmark",
    "Trøndelag"
  )
  
  Encoding(counties) <- "UTF-8"
  return(counties)
}
