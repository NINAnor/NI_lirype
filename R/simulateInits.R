#' Simulate values for initializing integrated model
#'
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#' @param R_perF logical. If TRUE, treats recruitment rate as juvenile per adult female.
#' If FALSE, treats recruitment rate as juvenile per adult (sum of both sexes).
#' @param shareRE logical. If TRUE, temporal random effects are shared across locations.
#' @param survVarT logical. If TRUE, survival is simulated including annual variation.
#' @param fitRodentCov logical. If TRUE, initial values are generated for rodent 
#' covariate (effect) and covariate is included in data simulation.
#' 
#' @return A list containing one complete set of initial values for the model.
#' @export
#'
#' @examples

simulateInits <- function(nim.data, nim.constants, R_perF, shareRE, survVarT, fitRodentCov){
  
  # Limits and constants #
  #----------------------#
  
  N_areas <- nim.constants$N_areas
  N_ageC <- nim.constants$N_ageC
  N_years <- nim.constants$N_years
  
  if(N_areas == 1){
    N_sites <- nim.constants$N_sites[1]
  }else{
    N_sites <- nim.constants$N_sites
  }
  
  L <- nim.data$L
  W <- nim.constants$W
  pi <- 3.141593
  A <- nim.data$A
  
  
  # Missing covariate values #
  #--------------------------#
  
  if(fitRodentCov) {
    if(!is.null(nim.data$RodentOcc)){
      RodentOcc <- nim.data$RodentOcc
    }else{
      RodentOcc <- matrix(0, nrow = N_areas, ncol = N_years)
    }
    
    if(NA %in% RodentOcc){
      RodentOcc[which(is.na(RodentOcc))] <- runif(length(which(is.na(RodentOcc))), 0, 1)
    }
    
    Inits_RodentOcc <- RodentOcc
    Inits_RodentOcc[which(!is.na(nim.data$RodentOcc))] <- NA
  }

  
  # Vital rates #
  #-------------#
  
  ## Area-specific survival parameters
  h.Mu.S <- runif(1, 0.3, 0.5)
  h.sigma.S <- runif(1, 0, 0.1)
  Mu.S1 <- runif(1, 0.6, 0.7)
  
  mu.S <-
    EnvStats::rnormTrunc(N_areas, qlogis(h.Mu.S), sd = h.sigma.S, max = qlogis(Mu.S1))
  Mu.S <- rep(NA, N_areas)
  S <-  matrix(NA, nrow = N_areas, ncol = N_years)
  
  epsT.S <- matrix(0, nrow = N_areas, ncol = N_years)
  #epsT.S <- matrix(rnorm(N_areas*N_years, 0, sigmaS.R), nrow = N_areas, ncol = N_years)
  
  for (x in 1:N_areas) {
    Mu.S[x] <- plogis(mu.S[x])
    S[x, 1:N_years] <-
      plogis(qlogis(Mu.S[x]) + epsT.S[x, 1:N_years])
  }
  
  if (survVarT) {
    sigmaT.S <- runif(1, 0, 0.05)
    S1 <- S2 <- rep(NA, N_years)
    epsT.S1.prop <- runif(1, 0.3, 0.8)
    
    S1[1:N_years] <-
      plogis(qlogis(Mu.S1) + epsT.S1.prop * epsT.S[nim.constants$SurvAreaIdx, 1:N_years])
    S2[1:N_years] <-
      S[nim.constants$SurvAreaIdx, 1:N_years] / S1[1:N_years]
  }
  
  
  ## Area-specific reproductive parameters
  h.Mu.R  <- runif(1, 1, 4)
  h.sigma.R <- runif(1, 0, 0.05)
  
  h.Mu.betaR.R <- runif(1, 0, 0.1)
  h.sigma.betaR.R <- runif(1, 0, 0.05)
  
  Mu.R <- rlnorm(N_areas, meanlog = log(h.Mu.R), sdlog =  h.sigma.R)
  
  if(fitRodentCov){
    betaR.R <- rnorm(N_areas, mean = h.Mu.betaR.R, sd = h.sigma.betaR.R)
  }else{
    betaR.R <- rep(0, N_areas)
  }
    
  sigmaT.R <- runif(1, 0, 2)
  
  R_year <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  if(shareRE){
    #epsT.R <- rep(0, N_years)
    epsT.R <- rnorm(N_years, 0, sigmaT.R)
    
    for(x in 1:N_areas){
      if(fitRodentCov) {
        R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R[x]*RodentOcc[x, 1:N_years] + epsT.R[1:N_years])
      } else {
        R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R[x] + epsT.R[1:N_years])
      }
    }
  }else{
    #epsT.R <- matrix(0, nrow = N_areas, ncol = N_years)
    epsT.R <- matrix(rnorm(N_areas*N_years, 0, sigmaT.R), nrow = N_areas, ncol = N_years)
    
    for(x in 1:N_areas){
      if(fitRodentCov) {
        R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R[x]*RodentOcc[x, 1:N_years] + epsT.R[x, 1:N_years])
      } else {
        R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R[x] + epsT.R[x, 1:N_years])
      }
    }
  }

  
  # Detection parameters #
  #----------------------#
  
  ## Area-specific detection parameters
  h.mu.dd <- runif(1, 2, 5)
  h.sigma.dd <- runif(1, 0, 1)
  
  #mu.dd <- rnorm(N_areas, h.mu.dd, sd = h.sigma.dd)
  mu.dd <- rep(h.mu.dd, N_areas)
  sigmaT.dd <- runif(1, 1, 10)
  
  sigma <- esw <- p <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  
  if(shareRE){
    epsT.dd <- rep(0, N_years)
    #epsT.dd <- rnorm(N_years, 0, sd = sigmaT.dd)
    
    for(x in 1:N_areas){
      sigma[x, 1:N_years] <- exp(mu.dd[x] + epsT.dd[1:N_years])
    }
    
  }else{
    epsT.dd <- matrix(0, nrow = N_areas, ncol = N_years)
    #epsT.dd <- matrix(rnorm(N_areas*N_years, 0, sigmaT.dd), nrow = N_areas, ncol = N_years)
    
    for(x in 1:N_areas){
      sigma[x, 1:N_years] <- exp(mu.dd[x] + epsT.dd[x, 1:N_years])
    }
  }

  
  for(x in 1:N_areas){
    esw[x, 1:N_years] <- sqrt(pi * sigma[x, 1:N_years]^2 / 2) 
    p[x, 1:N_years] <- min(esw[x, 1:N_years], W) / W
  }
  
  
  # Population model #
  #------------------#
  
  ## Initial densities / population sizes
  Mu.D1 <- rep(NA, N_areas)
  sigma.D <- runif(N_areas, 0.05, 2)
  ratio.JA1 <- runif(N_areas, 0.2, 0.6)
  
  N_exp <- Density <-
    array(0, dim = c(N_areas, N_ageC, max(N_sites), N_years))
  
  for (x in 1:N_areas) {
    D_x_sum <-
      apply(nim.data$N_a_line_year[x, , , ], c(2, 3), sum) / (L[x, , ] * W *
                                                                2)
    D_data <- D_x_sum[which(!is.na(D_x_sum) & D_x_sum > 0)]
    Mu.D1[x] <-
      runif(1, quantile(D_data, 0.25), quantile(D_data, 0.75))
    
    for (j in 1:N_sites[x]) {
      Density[x, 1, j, 1] <- Mu.D1[x] * ratio.JA1[x]
      Density[x, 2, j, 1] <- Mu.D1[x] * (1 - ratio.JA1[x])
      
      N_exp[x, 1, j, 1] <-
        extraDistr::rtpois(1, lambda = Density[x, 1, j, 1] * L[x, j, 1] * W * 2, a = 1)
      N_exp[x, 2, j, 1] <-
        extraDistr::rtpois(1, lambda = Density[x, 2, j, 1] * L[x, j, 1] * W * 2, a = 1)
    }
  }

  ## Population projection over time
  for (x in 1:N_areas) {
    for (j in 1:N_sites[x]) {
      for (t in 2:N_years) {
        Density[x, 2, j, t] <-
          sum(Density[x, 1:N_ageC, j, t - 1]) * S[x, t - 1] # Adults
        
        if (R_perF) {
          Density[x, 1, j, t] <-
            (Density[x, 2, j, t] / 2) * R_year[x, t] # Juveniles
        } else{
          Density[x, 1, j, t] <- Density[x, 2, j, t] * R_year[x, t] # Juveniles
        }
        
        N_exp[x, 1:N_ageC, j, t] <-
          Density[x, 1:N_ageC, j, t] * L[x, j, t] * W * 2
      }
    }
  }
  
  ## Area-specific population size and density
  N_tot_exp <- matrix(NA, nrow = N_areas, ncol = N_years)
  
  for (x in 1:N_areas) {
    for (t in 1:N_years) {
      N_tot_exp[x, t] <-
        sum(N_exp[x, 1, 1:N_sites[x], t] + N_exp[x, 2, 1:N_sites[x], t])    ## Summing up expected number of birds in covered area;
    }
  }
  
  # Assembly #
  #----------#
  
  InitVals <- list(
    #b = runif(1, 1, 50), 
    
    Mu.D1 = Mu.D1, 
    sigma.D = sigma.D,
    eps.D1 = matrix(0, nrow = nim.constants$N_areas, ncol = max(N_sites)),
    ratio.JA1 = ratio.JA1,
    
    Mu.R = Mu.R,
    h.Mu.betaR.R = h.Mu.betaR.R, h.sigma.betaR.R = h.sigma.betaR.R,
    h.Mu.R = h.Mu.R, h.sigma.R = h.sigma.R,
    sigmaT.R = sigmaT.R,
    epsT.R = epsT.R, 
    R_year = R_year,
    
    mu.dd = mu.dd,
    h.mu.dd = h.mu.dd, h.sigma.dd = h.sigma.dd,
    sigmaT.dd = sigmaT.dd,
    epsT.dd = epsT.dd,
    sigma = sigma, sigma2 = sigma^2,
    esw = esw,
    p = p,
    
    h.Mu.S = h.Mu.S,
    h.sigma.S = h.sigma.S,
    Mu.S1 = Mu.S1,
    mu.S = mu.S,
    Mu.S = Mu.S,
    epsT.S = epsT.S,
    S = S,
    
    Density = Density,
    N_exp = N_exp,
    N_tot_exp = N_tot_exp
  )
  
  if(survVarT){
    InitVals$sigmaT.S <- sigmaT.S
    InitVals$epsT.S1.prop <- epsT.S1.prop
    InitVals$S1 <- S1
    InitVals$S2 <- S2
  }
  
  if(fitRodentCov){
    InitVals$h.Mu.betaR.R <- h.Mu.betaR.R
    InitVals$h.sigma.betaR.R <- h.sigma.betaR.R
    InitVals$betaR.R <- betaR.R
    InitVals$RodentOcc <- Inits_RodentOcc
  }
  
  return(InitVals)
}
