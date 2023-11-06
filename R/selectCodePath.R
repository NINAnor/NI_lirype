#' Select appropriate code file based on toggle values
#'
#' @param shareRE logical. If TRUE, temporal random effects are shared across
#' locations.
#' @param survVarT logical. If TRUE, survival is modeled including among-year
#' variation.
#' @param addDummyDim logical. If TRUE (default) additional location/area
#' dimension is used.
#'
#' @return
#' @export
#'
#' @examples

selectCodePath <- function(shareRE, survVarT, addDummyDim = TRUE){
  
  if(addDummyDim){
    
    if(survVarT){
      code.path <- "NIMBLE code/rypeIDSM_multiArea_dHN_sepRE_survT.R"
    }else{
      if(shareRE){
        code.path <- "NIMBLE code/rypeIDSM_multiArea_dHN.R"
      }else{
        code.path <- "NIMBLE code/rypeIDSM_multiArea_dHN_sepRE.R"
      }
    }
    
  }else{
    
    if(survVarT){
      code.path <- "NIMBLE code/rypeIDSM_dHN_survT.R"
    }else{
      code.path <- "NIMBLE code/rypeIDSM_dHN.R"
    }
  }
  
  return(code.path)
}