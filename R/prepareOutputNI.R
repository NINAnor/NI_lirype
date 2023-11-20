#' Prepare relevant output data from NIMBLE model for subsequent implementation
#' in the NI database
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param N_areas integer. Number of areas included in the analyses.
#' @param area_names character vector containing area/location names.
#' @param N_sites matrix containing the number of sites per location/area/county.
#' @param min_years integer vector. Indices of first year with available data
#' for each location/area/county.
#' @param max_years integer vector. Indices of last year with available data for
#' each location/area/county.
#' @param minYear integer. The first year considered in analyses.
#' @param maxYear integer. The last year considered in analyses.
#'
#' @return a list of data frames with density estimates:
#'   total = average density estimates per area for the whole data period
#'   total.last10yr = average density estimates per area for the last 10 years of data
#'   annual.area = density estimates per area per year
#'   annual.area.adult = density estimates of adult ptarmigan per area per year
#' @export
#'
#' @examples

prepareOutputNI <- function(mcmc.out,
                            N_areas,
                            area_names,
                            N_sites,
                            min_years,
                            max_years,
                            minYear,
                            maxYear) {
  
  ## Convert posterior samples to matrix format
  out.mat <- as.matrix(mcmc.out)
  
  ## Summarize posteriors for population density
  popDens.areayear <- data.frame() # density per area and year
  popDens.ad.areayear <- data.frame() #adult density per area and year
  popDens.tot <- data.frame() #total density for the whole time period
  popDens.10yr <- data.frame() #total density for the last 10 years
  
  ## Annual densities per area
  for (i in 1:N_areas) {
    # Prepare matrix for temporary storage of results
    popDens.sum <-  data.frame()
    popDens.ad.sum <- data.frame()
    
    # Determine area-specific year range
    area_yearIdxs <- (1:(maxYear - minYear + 1))
    area_years <- area_yearIdxs + (minYear - 1)
    
    for (t in 1:length(area_years)) {
      # Summarize annual average population densities
      popDens_juv <- out.mat[, paste0("Density[",  i, ", 1, ", 1:N_sites[i], ", ", area_yearIdxs[t], "]")]
      popDens_ad <- out.mat[, paste0("Density[",  i, ", 2, ", 1:N_sites[i], ", ", area_yearIdxs[t], "]")]
      
      if (N_sites[i] > 1) {
        popDens_mean <- rowMeans(popDens_juv + popDens_ad)
      } else{
        popDens_mean <- popDens_juv + popDens_ad
      }
      
      # Merge average density back into posterior samples (for later calculation of total densities)
      out.mat <- cbind(out.mat, popDens_mean)
      colnames(out.mat)[dim(out.mat)[2]] <- paste0("meanDens[", i, ", ", area_yearIdxs[t], "]")
      
      # Mean densities per year
      popDens_add <- data.frame(
        Area = area_names[i],
        Year = area_years[t],
        Median = median(popDens_mean),
        lCI = unname(quantile(popDens_mean, probs = 0.25)), #25% lower
        uCI = unname(quantile(popDens_mean, probs = 0.75)), #75% upper
        Mean = mean(popDens_mean),
        SD = sd(popDens_mean)
      )
      popDens.sum <- rbind(popDens.sum, popDens_add)
      
      # Mean adult densities per year
      popDens_add.ad <- data.frame(
        Area = area_names[i],
        Year = area_years[t],
        Median = median(popDens_ad),
        lCI = unname(quantile(popDens_ad, probs = 0.25)), #25% lower
        uCI = unname(quantile(popDens_ad, probs = 0.75)), #75% upper
        Mean = mean(popDens_ad),
        SD = sd(popDens_ad)
      )
      popDens.ad.sum <- rbind(popDens.ad.sum, popDens_add.ad)
    }
    
    popDens.areayear <- rbind(popDens.areayear, popDens.sum)
    popDens.ad.areayear <- rbind(popDens.ad.areayear, popDens.ad.sum)
  }
  
  ## Calculate average area-specific densities over different time periods
  for (i in 1:N_areas) {
    # Entire time period - average total density
    dens_mean_tot <- rowMeans(out.mat[, paste0("meanDens[", i, ", ", area_yearIdxs[1:(length(area_years))], "]")])
    out.mat <- cbind(out.mat, dens_mean_tot)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("densAvg_tot[", i, "]")
    
    # last 10 years - average total density
    dens_mean_10yr <- rowMeans(out.mat[, paste0("meanDens[", i, ", ", area_yearIdxs[(length(area_years) - 10):length(area_years)], "]")])
    out.mat <- cbind(out.mat, dens_mean_10yr)
    colnames(out.mat)[dim(out.mat)[2]] <- paste0("densAvg_10yr[", i, "]")
  }
  
  ## Summarize average population densities
  for (i in 1:N_areas) {
    Dens_names <-
      c(paste0("densAvg_tot[",  i, "]"),
        paste0("densAvg_10yr[",  i, "]"))
    Dens_tot_add <- data.frame(
      Area = area_names[i],
      SummaryPeriod = paste0(minYear, "-", maxYear),
      Median = median(out.mat[, Dens_names[1]]),
      lCI = unname(quantile(out.mat[, Dens_names[1]], probs = 0.25)), #25% lower
      uCI = unname(quantile(out.mat[, Dens_names[1]], probs = 0.75)), #75% upper
      Mean = mean(out.mat[, Dens_names[1]]),
      SD = sd(out.mat[, Dens_names[1]])
    )
    Dens_10yr_add <- data.frame(
      Area = area_names[i],
      SummaryPeriod = paste0(maxYear - 10, "-", maxYear),
      Median = median(out.mat[, Dens_names[2]]),
      lCI = unname(quantile(out.mat[, Dens_names[2]], probs = 0.25)), #25% lower
      uCI = unname(quantile(out.mat[, Dens_names[2]], probs = 0.75)), #75% upper
      Mean = mean(out.mat[, Dens_names[2]]),
      SD = sd(out.mat[, Dens_names[2]])
    )
    popDens.tot <- rbind(popDens.tot, Dens_tot_add)
    popDens.10yr <- rbind(popDens.10yr, Dens_10yr_add)
  }
  
  return(
    list(
      total = popDens.tot,
      total.last10yr = popDens.10yr,
      annual.area = popDens.areayear,
      annual.area.adult = popDens.ad.areayear
    )
  )
}