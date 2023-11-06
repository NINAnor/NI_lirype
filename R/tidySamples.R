#' Tidy up posterior samples from multi-area model
#' 
#' Rescales density measures from square metres (as implemented in the model) to 
#' square kilometres (commonly used for reporting). 
#' It also removes all redundant nodes. These are nodes within the population 
#' size and density arrays with area-site-year combinations that do not appear 
#' in the data and are not estimated in the model. The function identifies these
#' nodes by the presence of NAs / only zeros in the posterior samples. 
#'
#' @param IDSM.out an MCMC list containing posterior samples for one or more
#' chains from a fitted multi-area model.
#' @param save logical. If TRUE (default), saves the output in an .RDS file.
#'
#' @return a condensed version of the MCMC list `IDSM.out` without redundant nodes.
#' @export
#'
#' @examples
#' 
tidySamples <- function(IDSM.out,
                        save = TRUE) {
  
  ## Extract number of chains
  n_chains <- length(IDSM.out)
  
  ## Identify density nodes
  Dnode_idx <- which(stringr::str_detect(colnames(IDSM.out[[1]]), pattern = "Density"))
  
  ## Identify unnecessary NA nodes
  NAnode_idx <- unname(which(is.na(colSums(IDSM.out[[1]]))))
  
  ## Identify 0 abundance nodes
  ZEROnode_idx <- unname(which(colSums(IDSM.out[[1]]) == 0))
  
  ## Combine into removal nodes
  REMOVEnode_idx <- c(NAnode_idx, ZEROnode_idx)
  
  for(i in 1:n_chains){
    
    ## Scale density measures (m^2 --> km^2)
    IDSM.out[[i]][, Dnode_idx] <- IDSM.out[[i]][, Dnode_idx]*1000^2
    
    ## Remove unnecessary NA & 0 nodes
    if(length(REMOVEnode_idx) != 0){
      IDSM.out[[i]] <- IDSM.out[[i]][, -REMOVEnode_idx]
    }
  }
  
  ## Optional: save tidied samples
  if(save){
    saveRDS(IDSM.out, file = "rypeIDSM_dHN_multiArea_realData_allAreas_tidy.rds")
  }
  
  ## Return tidied samples
  return(IDSM.out)
}