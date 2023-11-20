#' Plot MCMC traces and posterior densities to pdf
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param fitRodentCov logical. If TRUE, rodent covariate on reproduction is included.
#' @param VitalRates logical. If TRUE (default), plots traces and posterior densities for vital rate parameters.
#' @param DetectParams logical. If TRUE (default), plots traces and posterior densities for detection parameters.
#' @param PopSizes logical. If TRUE (default), plots traces and posterior densities for population sizes.
#' @param Densities logical. If TRUE (default), plots traces and posterior densities for population densities.
#'
#' @return a vector of pdf plot names. The plots can be found in Plots/MCMCTraces.
#' @export
#'
#' @examples
#' 
plotMCMCTraces <- function(mcmc.out, fitRodentCov, VitalRates = TRUE, DetectParams = TRUE, PopSizes = TRUE, Densities = TRUE){
  
  ## Make parameter lists
  mVR_params <- c(
    "Mu.R",
    "h.Mu.R",
    "h.sigma.R",
    "sigmaT.R",
    "Mu.D1",
    "sigma.D",
    "ratio.JA1",
    "Mu.S",
    "Mu.S1",
    "h.Mu.S",
    "h.sigma.S"
  )
  
  if (fitRodentCov) {
    mVR_params <-
      c(mVR_params, "betaR.R", "h.Mu.betaR.R", "h.sigma.betaR.R")
  }
  
  tVR_params <- c("R_year", "S")
  
  mDet_params <- c("mu.dd", "sigmaT.dd")
  
  tDet_params <- c("esw", "p", "sigma")
  
  N_params <- c("N_exp")
  
  Ntot_params <- c("N_tot_exp")
  
  Dens_params <- c("Density")
  
  
  ## Make plots and print to pdf
  plot.paths <- c()
  # Check if folder exists, if not create folder (NB! must use 'recursive' flag for first run)
  ifelse(!dir.exists("Plots/MCMCTraces"), dir.create("Plots/MCMCTraces", recursive = T), FALSE)

  # Vital rates
  if(VitalRates) {
    pdf("Plots/MCMCTraces/MCMCtrace_mVR.pdf",
        width = 8,
        height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = mVR_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
    
    pdf("Plots/MCMCTraces/MCMCtrace_tVR.pdf",
        width = 8,
        height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = tVR_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
    
    plot.paths <-
      c(
        plot.paths,
        "Plots/MCMCTraces/MCMCtrace_mVR.pdf",
        "Plots/MCMCTraces/MCMCtrace_tVR.pdf"
      )
  }
  
  # Detection parameters
  if(DetectParams) {
    pdf("Plots/MCMCTraces/MCMCtrace_mDet.pdf",
        width = 8,
        height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = mDet_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
    
    pdf("Plots/MCMCTraces/MCMCtrace_tDet.pdf",
        width = 8,
        height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = tDet_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
    
    plot.paths <-
      c(
        plot.paths,
        "Plots/MCMCTraces/MCMCtrace_mDet.pdf",
        "Plots/MCMCTraces/MCMCtrace_tDet.pdf"
      )
  }
  
  # Population sizes
  if(PopSizes) {
    pdf("Plots/MCMCTraces/MCMCtrace_Ntot.pdf",
        width = 8,
        height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = Ntot_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
    
    pdf("Plots/MCMCTraces/MCMCtrace_N.pdf",
        width = 8,
        height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = N_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
    
    plot.paths <-
      c(
        plot.paths,
        "Plots/MCMCTraces/MCMCtrace_Ntot.pdf",
        "Plots/MCMCTraces/MCMCtrace_N.pdf"
      )
  }
  
  # Population densities
  if(Densities) {
    pdf("Plots/MCMCTraces/MCMCtrace_Dens.pdf",
        width = 8,
        height = 6)
    MCMCvis::MCMCtrace(mcmc.out,
                       params = Dens_params,
                       ind = TRUE,
                       pdf = FALSE)
    dev.off()
    
    plot.paths <-
      c(plot.paths, "Plots/MCMCTraces/MCMCtrace_Dens.pdf")
  }
  
  return(plot.paths)
}
