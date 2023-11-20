rypeIDSM <- nimbleCode({
  
  ####################################################
  #### Distance sampling half normal detection function; 
  ## Data; 
  
  # y <- vector with distances to transect line; 
  # N_years <- number of years
  # N_obs <- number of observations
  # Year_obs <- vector with year for each observation
  
  
  # priors for distance model; 
  
  pi <- 3.141593
  
  #random year effect for distance sampling model; 
  ## Priors for hyper-parameters
  h.mu.dd ~ dunif(-10, 100)
  h.sigma.dd ~ dunif(0, 5)

  for(x in 1:N_areas){
    mu.dd[x] ~ dnorm(h.mu.dd, sd = h.sigma.dd)
  }
  
  sigmaT.dd ~ dunif(0, 20)
  
  for(t in 1:N_years){
    epsT.dd[t] ~ dnorm(0, sd = sigmaT.dd)
  }
  #TODO: Implement spatial correlation in temporal RE (epsT.dd[x, t])
  
  
  ########################################################
  for(x in 1:N_areas){
    for(t in 1:N_years){
      
      log(sigma[x, t]) <- mu.dd[x] + epsT.dd[t]
      #log(sigma[x, t]) <- mu.dd[x] + epsT.dd[x, t]
      
      sigma2[x, t] <- sigma[x, t] * sigma[x, t]
      
      # effective strip width
      esw[x, t] <- sqrt(pi * sigma2[x, t] / 2) 
      p[x, t] <- min(esw[x, t], W) / W
    }
  }
  
  ########################################################   
  for(x in 1:N_areas){
    for (i in 1:N_obs[x]){ 
      # LIKELIHOOD
      # using nimbleDistance::dHN
      y[x, i] ~ dHN(sigma = sigma[x, Year_obs[x, i]], Xmax = W, point = 0)
    }
  }

  ###################################################
  ## Random effects model for R (i.e. )
  ## Data: 
  
  # R_obs <- vector with number of chicks / obs [0 - 12]
  # N_Years <- number of years in time series 
  # N_obs <- number of observations 
  # Year_obs <- vector with years for each observation. 
  
  ## Priors; 

  h.Mu.R  ~ dunif(0, 15)
  h.sigma.R ~ dunif(0, 5)
  
  if(fitRodentCov){
    h.Mu.betaR.R ~ dunif(-5, 5)
    h.sigma.betaR.R ~ dunif(0, 5)
  }
  
  sigmaT.R ~ dunif(0, 5)
  
  for (t in 1:N_years){
    epsT.R[t] ~ dnorm(0, sd = sigmaT.R) # Temporal RE (shared across areas)
  }
  #TODO: Implement spatial correlation in temporal RE (epsT.R[x, t])
  

  for(x in 1:N_areas){
    
    Mu.R[x]  ~ dlnorm(meanlog = log(h.Mu.R), sdlog = h.sigma.R)
    
    if(fitRodentCov){
      betaR.R[x] ~ dnorm(mean = h.Mu.betaR.R, sd = h.sigma.betaR.R)
    }
    
    ## Constraints;
    if(fitRodentCov){
      R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + betaR.R[x]*RodentOcc[x, 1:N_years] + epsT.R[1:N_years])
    }else{
      R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + epsT.R[1:N_years])
      #R_year[x, 1:N_years] <- exp(log(Mu.R[x]) + epsT.R[x, 1:N_years])
    }
    
    ## Likelihood;
    for (i in 1:N_sumR_obs[x]){
      sumR_obs[x, i] ~ dpois(R_year[x, sumR_obs_year[x, i]]*sumAd_obs[x, i])
    }
  }
  
  
  ########################################################################
  ### MODEL FOR Density in year 1:
  ### Simple random effects model 
  # Data; 
  
  # N_sites <- number of sites
  # N_line_year <- number of birds pr. line 
  # L <- length of transect lines
  
  for(x in 1:N_areas){
    
    ## Priors; 
    
    for(j in 1:N_sites[x]){
      eps.D1[x, j] ~ dnorm(0, sd = sigma.D[x])
    }
    
    Mu.D1[x] ~ dunif(0, 10)
    sigma.D[x] ~ dunif(0, 20)
    
    ratio.JA1[x] ~ dunif(0, 1)
    
    ## State model
    for (j in 1:N_sites[x]){
      
      #for(a in 1:N_ageC){
      #  N_exp[a, j, 1] ~ dpois(Density[a, j, 1]*L[1, j, 1]*W*2)      ## Expected number of birds
      #}  
      
      N_exp[x, 1, j, 1] ~ dpois(Density[x, 1, j, 1]*L[x, j, 1]*W*2) 
      N_exp[x, 2, j, 1] ~ dpois(Density[x, 2, j, 1]*L[x, j, 1]*W*2) 
      
      Density[x, 1, j, 1] <- exp(log(Mu.D1[x]) + eps.D1[x, j])*ratio.JA1[x]             ## random effects model for spatial variation in density for year 1
      Density[x, 2, j, 1] <- exp(log(Mu.D1[x]) + eps.D1[x, j])*(1-ratio.JA1[x])
      
      ## Detection model year 1
      for(a in 1:N_ageC){
        N_a_line_year[x, a, j, 1] ~ dpois(p[x, 1]*N_exp[x, a, j, 1])
      }
      
      #N_line_year[x, j, 1] ~ dpois(p[x, 1]* sum(N_exp[x, 1:N_ageC, j, 1]))
    }
  }
  
  #####################################################
  ## Model for survival; 
  
  # Priors
  h.Mu.S ~ dunif(0, 1)
  h.sigma.S ~ dunif(0, 5)

  Mu.S1 ~ dunif(0, 1) # Season 1 survival for Lierne

  for(x in 1:N_areas){
    mu.S[x] ~ dnorm(logit(h.Mu.S), sd = h.sigma.S)
    #TODO: Implement spatial correlation in survival averages

    ## Constraints
    logit(Mu.S[x]) <- mu.S[x]
    S[x, 1:N_years] <- Mu.S[x]
    #TODO: Consider implementing spatially correlated random year variation

  }

  if (survVarT) {
    S1[1:N_years] <- Mu.S1
    S2[1:N_years] <- S[SurvAreaIdx, 1:N_years] / S1[1:N_years]
    
    ## Data likelihoods
    for (t in 1:N_years_RT) {
      Survs1[t, 2] ~ dbinom(S1[year_Survs[t]], Survs1[t, 1])
      Survs2[t, 2] ~ dbinom(S2[year_Survs[t]], Survs2[t, 1])
    }
  }
  
  
  #####################################################    
  ### Model for year 2 - n.years; 
  ### post-breeding census
  
  for (x in 1:N_areas) {
    for (j in 1:N_sites[x]) {
      for (t in 2:N_years) {
        ## Process model
        Density[x, 2, j, t] <-
          sum(Density[x, 1:N_ageC, j, t - 1]) * S[x, t - 1] # Adults
        
        if (R_perF) {
          Density[x, 1, j, t] <-
            (Density[x, 2, j, t] / 2) * R_year[x, t] # Juveniles
        } else {
          Density[x, 1, j, t] <- Density[x, 2, j, t] * R_year[x, t] # Juveniles
        }
        
        N_exp[x, 1:N_ageC, j, t] <-
          Density[x, 1:N_ageC, j, t] * L[x, j, t] * W * 2
        
        ## Detection model year 2 - T
        for (a in 1:N_ageC) {
          N_a_line_year[x, a, j, t] ~ dpois(p[x, t] * N_exp[x, a, j, t])
        }
        
        #N_line_year[x, j, t] ~ dpois(p[x, t]*sum(N_exp[x, 1:N_ageC, j, t]))
      }
    }
  }
  
  ####################################################
  ## Observation model
  ## P is estimated in distance sampling component - based on 
  ## distance to transect line data; 
  
  ####################################################
  ### Derived parameters; Nt and Dt
  for (x in 1:N_areas) {
    for (t in 1:N_years) {
      N_tot_exp[x, t] <- sum(N_exp[x, 1, 1:N_sites[x], t] + N_exp[x, 2, 1:N_sites[x], t])    ## Summing up expected number of birds in covered area;
      #D[x, t] <- N_tot_exp[x, t] / A[x, t]       ## Deriving density as N/A
    }
  }
  
  ####################################################
  ## Rodent covariate imputation model
  if(fitRodentCov){
    for(x in 1:N_areas){
      for (t in 1:N_years){
        RodentOcc[x, t] ~ dunif(0, 1)
      }
    }
  }
  
})
