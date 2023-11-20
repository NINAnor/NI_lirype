#' Plot time series of vital rates, detection probabilities, and population densities to pdf
#'
#' @param mcmc.out an mcmc list containing posterior samples from a model run.
#' @param N_areas integer. Number of areas included in analyses. 
#' @param area_names character vector containing area/location names.
#' @param N_sites matrix containing the number of sites per area/location.
#' @param min_years integer vector. Indices of first year with available data for each area/location. 
#' @param max_years integer vector. Indices of last year with available data for each area/location. 
#' @param minYear integer. The first year considered in analyses. 
#' @param maxYear integer. The last year considered in analyses. 
#' @param VitalRates logical. If TRUE (default), plots time series of vital rate parameters.
#' @param DetectParams logical. If TRUE (default), plots time series of detection parameters.
#' @param Densities logical. If TRUE (default), plots time series of average population densities.
#' @param survVarT logical. Global switch for survival estimates.
#' @param showDataWindow logical. If TRUE (default) highlights the data collection period in plots.
#'
#' @return a vector of pdf plot names. The plots can be found in Plots/TimeSeries.
#' @export
#'
#' @examples

plotTimeSeries <- function(mcmc.out, 
                           N_areas, area_names, N_sites, 
                           min_years, max_years, minYear, maxYear,
                           VitalRates = TRUE, DetectParams = TRUE, Densities = TRUE, survVarT,
                           showDataWindow = TRUE){
  
  ## Check that necessary packages are loaded and installed
  if(!require(ggplot2)) {
    print("Package 'ggplot2' not installed")
    break
  }
  
  ## Convert posterior samples to matrix format
  out.mat <- as.matrix(mcmc.out)
  
  ## Summarize posteriors for relevant parameters
  rRep <- data.frame()
  if (survVarT) {
    pSurv <- data.frame()
  }
  pDetect <- data.frame()
  popDens <- data.frame()
  
  for(i in 1:N_areas) {
    # Prepare matrices for temporary storage of results
    rRep.sum <- data.frame()
    if (survVarT) {
      pSurv.sum <- data.frame()
    }
    pDetect.sum <- data.frame()
    popDens.sum <-  data.frame()
    
    # Determine area-specific year range
    #area_yearIdxs <- (min_years[i]:max_years[i])
    area_yearIdxs <- (1:(maxYear - minYear + 1))
    area_years <- area_yearIdxs + (minYear - 1)
    
    for (t in 1:length(area_years)) {
      # Summarize annual reproductive rates
      rRep_name <- paste0("R_year[",  i, ", ", area_yearIdxs[t], "]")
      rRep_add <- data.frame(
        Area = area_names[i],
        Year = area_years[t],
        Median = median(out.mat[, rRep_name]),
        lCI = unname(quantile(out.mat[, rRep_name], probs = 0.025)),
        uCI = unname(quantile(out.mat[, rRep_name], probs = 0.975))
      )
      rRep.sum <- rbind(rRep.sum, rRep_add)
      
      # Summarize annual survival rates
      if (survVarT) {
        pSurv_name <- paste0("S[",  i, ", ", area_yearIdxs[t], "]")
        pSurv_add <- data.frame(
          Area = area_names[i],
          Year = area_years[t],
          Median = median(out.mat[, pSurv_name]),
          lCI = unname(quantile(out.mat[, pSurv_name], probs = 0.025)),
          uCI = unname(quantile(out.mat[, pSurv_name], probs = 0.975))
        )
        pSurv.sum <- rbind(pSurv.sum, pSurv_add)
      }
      
      # Summarize annual detection probabilities
      pDetect_name <- paste0("p[",  i, ", ", area_yearIdxs[t], "]")
      pDetect_add <- data.frame(
        Area = area_names[i],
        Year = area_years[t],
        Median = median(out.mat[, pDetect_name]),
        lCI = unname(quantile(out.mat[, pDetect_name], probs = 0.025)),
        uCI = unname(quantile(out.mat[, pDetect_name], probs = 0.975))
      )
      pDetect.sum <- rbind(pDetect.sum, pDetect_add)
      
      # Summarize annual average population densities
      popDens_juv <- out.mat[, paste0("Density[",  i, ", 1, ", 1:N_sites[i], ", ", area_yearIdxs[t], "]")]
      popDens_ad <- out.mat[, paste0("Density[",  i, ", 2, ", 1:N_sites[i], ", ", area_yearIdxs[t], "]")]
      #popDens_mean <- rowMeans(popDens_juv + popDens_ad) This gives error if number of transects within area is 1
      
      if (N_sites[i] > 1) {
        popDens_mean <- rowMeans(popDens_juv + popDens_ad)
      } else{
        popDens_mean <- popDens_juv + popDens_ad
      }
      
      popDens_add <- data.frame(
        Area = area_names[i],
        Year = area_years[t],
        Median = median(popDens_mean),
        lCI = unname(quantile(popDens_mean, probs = 0.025)),
        uCI = unname(quantile(popDens_mean, probs = 0.975))
      )
      popDens.sum <- rbind(popDens.sum, popDens_add)
    }
    
    rRep <- rbind(rRep, rRep.sum)
    if (survVarT) {
      pSurv <- rbind(pSurv, pSurv.sum)
    }
    pDetect <- rbind(pDetect, pDetect.sum)
    popDens <- rbind(popDens, popDens.sum)
  }
  

  ## Make plots and print to pdf
  plot.paths <- c()
  ifelse(
    !dir.exists("Plots/TimeSeries"),
    dir.create("Plots/TimeSeries", recursive = T),
    FALSE
  ) ## Check if folder exists, if not create folder
  
  # Reproductive rates
  if(VitalRates) {
    pdf("Plots/TimeSeries/TimeSeries_rRep.pdf",
        width = 8,
        height = 5)
    for (i in 1:N_areas) {
      p_rRep <-
        ggplot(subset(rRep, Area == area_names[i]), aes(x = Year)) +
        geom_line(aes(y = Median), color = "#67008A") +
        geom_ribbon(aes(ymin = lCI, ymax = uCI),
                    alpha = 0.5,
                    fill = "#67008A") +
        scale_x_continuous(breaks = c(minYear:maxYear),
                           limits = c(minYear, maxYear)) +
        ylim(min(rRep$lCI), max(rRep$uCI)) +
        ylab("Reproductive rate") +
        ggtitle(area_names[i]) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 0.75))
      
      if(showDataWindow) {
        p_rRep <- p_rRep  +
          geom_rect(
            xmin = min_years[i] + minYear - 1,
            xmax = max_years[i] + minYear - 1,
            ymin = min(rRep$lCI),
            ymax = max(rRep$uCI),
            alpha = 0.01,
            fill = "cornflowerblue"
          )
      }
      print(p_rRep)
      
    }
    dev.off()
    
    plot.paths <- c(plot.paths, "Plots/TimeSeries/TimeSeries_rRep.pdf")
  }
  
  # Survival probabilities
  if(survVarT) {
    pdf("Plots/TimeSeries/TimeSeries_pSurv.pdf",
        width = 8,
        height = 5)
    for (i in 1:N_areas) {
      p_pSurv <-
        ggplot(subset(pSurv, Area == area_names[i]), aes(x = Year)) +
        geom_line(aes(y = Median), color = "#9A42B8") +
        geom_ribbon(aes(ymin = lCI, ymax = uCI),
                    alpha = 0.5,
                    fill = "#9A42B8") +
        scale_x_continuous(#breaks = c(min_years[i]:max_years[i]) + minYear - 1,
          breaks = c(minYear:maxYear),
          limits = c(minYear, maxYear)) +
        ylim(min(pSurv$lCI), max(pSurv$uCI)) +
        ylab("Annual survival probability") +
        ggtitle(area_names[i]) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 0.75))
      
      if (showDataWindow) {
        p_pSurv <- p_pSurv +
          geom_rect(
            xmin = min_years[i] + minYear - 1,
            xmax = max_years[i] + minYear - 1,
            ymin = min(pSurv$lCI),
            ymax = max(pSurv$uCI),
            alpha = 0.01,
            fill = "cornflowerblue"
          )
      }
      print(p_pSurv)
      
    }
    dev.off()
    
    plot.paths <- c(plot.paths, "Plots/TimeSeries/TimeSeries_pSurv.pdf")
  }
  
  # Detection probabilities
  if(DetectParams) {
    pdf("Plots/TimeSeries/TimeSeries_pDetect.pdf",
        width = 8,
        height = 5)
    for (i in 1:N_areas) {
      p_pDetect <-
        ggplot(subset(pDetect, Area == area_names[i]), aes(x = Year))  +
        geom_line(aes(y = Median), color = "#856CEB") +
        geom_ribbon(aes(ymin = lCI, ymax = uCI),
                    alpha = 0.5,
                    fill = "#856CEB") +
        scale_x_continuous(#breaks = c(min_years[i]:max_years[i]) + minYear - 1,
          breaks = c(minYear:maxYear),
          limits = c(minYear, maxYear)) +
        ylim(min(pDetect$lCI), max(pDetect$uCI)) +
        ylab("Detection probability") +
        ggtitle(area_names[i]) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 0.75))
      
      if (showDataWindow) {
        p_pDetect <- p_pDetect +
          geom_rect(
            xmin = min_years[i] + minYear - 1,
            xmax = max_years[i] + minYear - 1,
            ymin = min(pDetect$lCI),
            ymax = max(pDetect$uCI),
            alpha = 0.01,
            fill = "cornflowerblue"
          )
      }
      print(p_pDetect)
      
    }
    dev.off()
    plot.paths <- c(plot.paths, "Plots/TimeSeries/TimeSeries_pDetect.pdf")
  }
  
  # Average population densities
  if(Densities) {
    pdf(
      "Plots/TimeSeries/TimeSeries_popDens1.pdf",
      width = 8,
      height = 5
    )
    for (i in 1:N_areas) {
      p_popDens1 <-
        ggplot(subset(popDens, Area == area_names[i]), aes(x = Year))  +
        geom_line(aes(y = Median), color = "#C2B391") +
        geom_ribbon(aes(ymin = lCI, ymax = uCI),
                    alpha = 0.5,
                    fill = "#C2B391") +
        scale_x_continuous(#breaks = c(min_years[i]:max_years[i]) + minYear - 1,
          breaks = c(minYear:maxYear),
          limits = c(minYear, maxYear)) +
        ylab(bquote("Average population density " (birds / km ^ 2))) +
        ggtitle(area_names[i]) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 0.75))
      
      if (showDataWindow) {
        p_popDens1 <- p_popDens1 +
          geom_rect(
            xmin = min_years[i] + minYear - 1,
            xmax = max_years[i] + minYear - 1,
            ymin = min(subset(popDens, Area == area_names[i])$lCI),
            ymax = max(subset(popDens, Area == area_names[i])$uCI),
            alpha = 0.01,
            fill = "cornflowerblue"
          )
      }
      print(p_popDens1)
      
    }
    dev.off()
    
    pdf(
      "Plots/TimeSeries/TimeSeries_popDens2.pdf",
      width = 8,
      height = 5
    )
    for (i in 1:N_areas) {
      p_popDens2 <-
        ggplot(subset(popDens, Area == area_names[i]), aes(x = Year))  +
        geom_line(aes(y = Median), color = "#C2B391") +
        geom_ribbon(aes(ymin = lCI, ymax = uCI),
                    alpha = 0.5,
                    fill = "#C2B391") +
        scale_x_continuous(#breaks = c(min_years[i]:max_years[i]) + minYear - 1,
          breaks = c(minYear:maxYear),
          limits = c(minYear, maxYear)) +
        ylim(min(popDens$lCI), max(popDens$uCI)) +
        ylab(bquote("Average population density " (birds / km ^ 2))) +
        ggtitle(area_names[i]) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 0.75))
      
      if (showDataWindow) {
        p_popDens2 <- p_popDens2 +
          geom_rect(
            xmin = min_years[i] + minYear - 1,
            xmax = max_years[i] + minYear - 1,
            ymin = min(popDens$lCI),
            ymax = max(popDens$uCI),
            alpha = 0.01,
            fill = "cornflowerblue"
          )
      }
      print(p_popDens2)
      
    }
    dev.off()
    
    plot.paths <-
      c(
        plot.paths,
        "Plots/TimeSeries/TimeSeries_popDens1.pdf",
        "Plots/TimeSeries/TimeSeries_popDens2.pdf"
      )
  }
  
  return(plot.paths)
}
