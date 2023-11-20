#' Download ptarmigan line transect data from Living Norway Data Portal
#'
#' @param datasets character (vector). The name(s) of the dataset to download.
#' Can be any combination of "Fjellstyrene", "Statskog", and "FeFo".
#' @param version numeric. Version(s) of the dataset(s) to download. 
#' @param save logical. If TRUE (default), saves downloaded .zip archive into 
#' "data" folder. Folder is created if it does not exist yet.
#'
#' @return A list of Darwin Core Archives (class DwCArchive, R6) containing line
#' transect data for willow ptarmigan (Lagopus lagopus) in Norway. 
#' @export
#'
#' @examples
#'      downloadLN(
#'        datasets = c("Fjellstyrene", "Statskog", "FeFo"),
#'        versions = c(1.7, 1.8, 1.12),
#'        save = FALSE
#'      )

downloadLN <- function(datasets, versions, save = TRUE){
  
  ## Define and check dataset names
  dataset_Names <- c("Fjellstyrene", "Statskog", "FeFo")
  if(any(!(datasets %in% dataset_Names))){
    stop("Invalid dataset name. Valid names include Fjellstyrene, Statskog, and FeFo.")
  }
  
  ## Define dataset keys
  dataset_Keys <- c("b49a2978-0e30-4748-a99f-9301d17ae119", # Fjellstyrene
                    "6a948a1c-7e23-4d99-b1c1-ec578d0d3159", # Statskog
                    "c47f13c1-7427-45a0-9f12-237aad351040") # FeFo
  
  ## Download specific data versions from Living Norway
  Rype_arkiv <- list()
  for(i in 1:length(datasets)){
    Rype_arkiv[[i]] <- LivingNorwayR::getLNportalData(datasetKey = dataset_Keys[which(dataset_Names == datasets[i])], version = versions[i])
  }
  names(Rype_arkiv) <- datasets 
  
  ## Optional: Save the archive file to ...data/
  if(save){
    dir.create("data", showWarnings = FALSE)
    for(i in 1:length(datasets)){
      Rype_arkiv[[i]]$exportAsDwCArchive(file.path(paste0("data/", datasets[i], "_arkiv.zip")))
    }
  }
  
  ## Return the archive
  return(Rype_arkiv)
}