#' Assign each transect line and observation to county by using the geographical
#' data available (WKT) for each transect line
#'
#' @param df a list of 2 tibbles with the wrangled line transect and observation
#' data
#'
#' @return a list of 2 tibbles with the assigned county name for transects and
#' observations
#' @export
#'
#' @examples
#'      assignCounty(LT_data,
#'                   counties = counties,
#'                   counties2020 = counties2020)

assignCounty <- function(df = LT_data, counties = NULL, counties2020) {
  #load required libraries
  if(!require(sf)) { print("Package 'sf' not installed"); break; }
  if(!require(dplyr)) { print("Package 'dplyr' not installed"); break; }
  
  df_trans <- df$d_trans
  df_obs <- df$d_obs
  
  #convert WKT to sf objects
  wgs4326 <- st_as_sfc(df_trans$footprintWKT, crs = 4326)

  #import county map
  if (!counties2020) {
    countyMap <- st_read("maps/counties.shp", quiet = T)
  } else {
    countyMap <- st_read("maps/counties2020.shp", quiet = T)
  }
  Encoding(countyMap$name) <- "latin1"
  countyMap <- countyMap[order(countyMap$OBJECTID),]
  rownames(countyMap) <- NULL
  
  #convert map attributes to data frame
  county.names <- as.data.frame(countyMap)
  county.names <- county.names %>% dplyr::select(OBJECTID, name)
  county.names$cname <- as.integer(rownames(county.names))

  #extract the county for each line transect
  countyIntersect <- st_intersects(wgs4326, countyMap, sparse = T)
  countyIntersect <- sapply(countyIntersect, "[[", 1) #retain only one county per transect line
  countyIntersect <- as.data.frame(countyIntersect)
  countyIntersect <- dplyr::left_join(countyIntersect, county.names, by = c("countyIntersect" = "cname"))
  colnames(countyIntersect)[3] <- "county"
  
  #assign county to transect lines
  df_trans$county <- countyIntersect$county
  
  #assign county to observations
  tmp.trans <- df_trans %>% dplyr::select(eventID, county)
  df_obs <- left_join(df_obs, tmp.trans, by = c("parentEventID" = "eventID"))
  
  #filter selected counties
  if(!is.null(counties)) {
    df_trans <- df_trans %>% dplyr::filter(county %in% counties)
    df_obs <- df_obs %>% dplyr::filter(county %in% counties)
  }
  
  return(list(d_trans = df_trans, d_obs = df_obs))
}
