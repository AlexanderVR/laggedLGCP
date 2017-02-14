#' Create grid and bin coalescent events
#' 
#' @param coal_data A list containing elements 'coal_times', 'samp_times', and 'n_sampled'.
#' @param n_bins Number of bins to use.
#' @param covariate If specified, model lag between eff. pop. size and covariate,
#' instead of eff. pop. size and sampling rate.
#' @param time_unit Unit of measurement for the timestamps.
#' 
#' @export
bin_coalescent_data <- function(coal_data, n_bins=100, covariate=NULL, time_unit='Years') {
  
  # Make sure event data is a list of numeric vectors
  stopifnot(all(vapply(coal_data, function(x) is.numeric(x) & is.vector(x), logical(1))))
  
  # Must contain coal_times, samp_times, and n_sampled
  stopifnot(all(sort(names(events)) == c('coal_times', 'n_sampled', 'samp_times')))
  
  st <- NULL
  fin <- NULL
  
  # if have covariate extend the grid to accomodate
  if (!is.null(covariate)) {
    st <- min(c(covariate, 0))
    fin <- max(c(covariate, coal_data$samp_times, coal_data$coal_times))
  }
  
  # Bin the event data
  pars <- bin_coal(coal_data, n_bins, st=st, fin=fin)
  
  # replace sampling times with covariate if specified
  if (!is.null(covariate)) {
    pars_cov <- bin_poisson(covariate, grid=pars$breaks)
    pars$counts <- cbind(pars$counts[, 1], pars_cov$counts)
    pars$offset <- cbind(pars$offset[, 1], pars_cov$offset) 
  }
  
  pars$time_unit <- time_unit
  
  # Indicate that the data is from a coalescent process
  pars$coalescent <- 1
  
  class(pars) <- "BinnedData"
  return(pars)
}


#' Create grid and bin events from multiple Poisson processes.
#' 
#' @param events A list of numeric vectors containing the timestamps.
#' @param n_bins Number of bins to use.
#' @param time_unit Unit of measurement for the timestamps.
#' @inheritParams 
#' @param fin End point of grid.
#' 
#' @export
bin_poisson_data <- function(events, n_bins=100, time_unit='Years', st=NULL, fin=NULL) {
  
  # Make sure event data is a list of numeric vectors
  stopifnot(all(vapply(events, function(x) is.numeric(x) & is.vector(x), logical(1))))
  
  # Bin the event data
  pars <- bin_poisson(events, n_bins, st=st, fin=fin)
  
  pars$time_unit <- time_unit
  
  # Indicate that the data is from Poisson process, not coalescent
  pars$coalescent <- 0
  
  class(pars) <- "BinnedData"
  return(pars)
}