#' Makes a list of (parent) eventIDs of duplicate transects that should be removed from data
#' 
#'Applies to the following locationIDs:
#'
#' locationID 18FD1A77-63FA-4DC0-AF40-E3934F68BC54, year 2016
#' --> Only one parentEventID has observations (BE03210E-1540-4106-A51B-D00EAC6C3F3D)
#' --> Drop other two: 3EAD98CE-C7C9-472D-9C79-48E735C4525D, 54F312E7-9E21-44C5-ADE3-90F47EC95F62
#'
#' locationID 34B7CD2A-1843-4234-8EEB-BBA9A9F460A7, year 2018
#' --> Perfect duplicate. Remove one. 
#' --> parentEventIDs: 10201BE7-2782-4062-B18D-0C7E17273821, B82C1256-D359-4075-B65B-4304CB086A4B
#'
#' locationID 3F324C69-3E3B-4502-BE5A-74048DF1C304, year 2016
#' --> Only one parentEventID has observations (480CE0F8-6311-4C42-9F1D-6CD9EA1FDD7D)
#' --> Drop other: 443FF8B8-071C-4AEB-B834-3C43748D194F
#'
#' locationID B1247B3C-0277-49DC-8714-4759A81E6E9B, year 2018
#' --> Perfect duplicate. Remove one. 
#' --> parentEventIDs: 37875623-B4AD-4651-A66C-2A2221AC6933, D7F9B1FA-75F9-43A7-9AD4-9479BB9B7C34
#'
#' locationID B6167D67-4448-4402-A26B-60303BE6CEF9, year 2018
#' --> Duplicate, but no observations. Remove one.  
#' --> parentEventIDs: 33B1EEEE-0130-48F5-9973-7FB196BB8A98, 20C98FA8-AB29-43D9-BC18-90973367B080
#'
#' locationID DA651143-0B57-44E1-8266-C9CC169391F2, year 2014
#' --> Duplicate, but no observations. Remove one.
#' --> parentEventIDs: 1B79AFEC-A512-4407-AC18-C77A890B4088, 05E5C283-74AB-4192-8E5F-E57370DDD527
#'
#' @return
#' @export
#'
#' @examples
#'     listDuplTransects()

listDuplTransects <- function(){

  ## List of (parent)EventIDs that are duplicates and should be dropped from the data
  duplTransects <- c(
   "3EAD98CE-C7C9-472D-9C79-48E735C4525D",
   "54F312E7-9E21-44C5-ADE3-90F47EC95F62",
   "B82C1256-D359-4075-B65B-4304CB086A4B",
   "443FF8B8-071C-4AEB-B834-3C43748D194F",
   "D7F9B1FA-75F9-43A7-9AD4-9479BB9B7C34",
   "20C98FA8-AB29-43D9-BC18-90973367B080",
   "05E5C283-74AB-4192-8E5F-E57370DDD527"
  )
  
  return(duplTransects)  
}
