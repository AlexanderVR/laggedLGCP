#' Generate a realization of log Gaussian Cox process, with each marginal latent function a shifted copy of the last.
#' 
#' @param n_series Number of correlated series.
#' @param lag The 'true' lag parameter, indicating how much each series is shifted to the right relative to the first series. Must be a vector of length = (n_series-1).
#' @param sim_res Resolution (bin-width) used to generate the time stamps.
#' @param type Which type of latent function to use. Either 'cyclic' or 'boombust'.
#' @param scaling Constant scaling factor for the overall Poisson rate.
#'
#' @export
sim_lag_poisson <- function(n_series = 2, lag = 0.5, sim_res=0.01, type='cyclic', scaling=.1) {
  # Generate realization of lagged poisson point processes.
  
  if (length(lag) != (n_series - 1))
    stop("length(lag) != n_series")
  if (type == 'cyclic')
    fun <- function(t) phylodyn::cyclic_traj(t)
  else if (type == 'boombust')
    fun <- function(t) phylodyn::boombust_traj(t, bust=15)
  else stop("unrecognized trajectory type")
  
  # extend lag vector to include first series
  lag = c(0, lag)

  # Set up grid for simulating the Poisson processes
  t <- seq(-40, 40, by=sim_res)
  dt <- t[2] - t[1]
  rates <- fun(t) * dt * scaling
  st <- -20
  fin <- 20
  events <- list()
  x <- c()
  for (k in 1:n_series) {
    y <- rpois(length(rates), rates)
    all_events <- t[y>0] + lag[k]
    events[[k]] <- all_events[(all_events > st) & (all_events < fin)]
    x <- cbind(x, rates[(t > (st - lag[k])) & (t < (fin - lag[k]))])
    
  }
  return(list(events=events, lag=lag, sim_res=sim_res, type=type, rates=x, t=t[t>st & t<fin]))
}

#' Bin a list of timestamps
#' 
#' @export
#' @param times List of vectors of timestamps.
#' @param n_bins Number of equally-spaced bins to use.
#' @param grid Existing grid to use.
#' @param st Assumed start time of data collection. Default = min(times) - 0.1 * bin_width
#' @param fin Assumed end time of data collection. Default = max(times) + 0.1 * bin_width
#'
bin_poisson <- function(times, n_bins=100, grid=NULL, st=NULL, fin=NULL) {
  # Input: list of vectors, one per timepoint
  if (is.numeric(times)) {
    times <- list(y=times)
  }
  n_series <- length(times)
  
  # Find first and last point
  if (is.null(grid)) {
    grid <- make_grid(times, n_bins, st=st, fin=fin)
  } else {
    n_bins <- length(grid) - 1
  }
  dt <- grid[2] - grid[1]
  offset <- matrix(dt, n_bins, n_series)
  counts <- matrix(0, n_bins, n_series)
  for (k in 1:n_series) {
    y <- times[[k]]
    counts[,k] <- hist(y, breaks=grid, plot=FALSE)$counts
  }
  return(list(counts=counts, breaks=grid, offset=offset, delta=dt, N=n_bins, 
              N_vis = n_bins * n_series, n_series=n_series))
}

make_grid <- function(times, n_bins, st, fin) {
  if (is.null(st))
    st <- min(vapply(times, FUN=min, FUN.VALUE=numeric(1)))
  if (is.null(fin))
    fin <- max(vapply(times, FUN=max, FUN.VALUE=numeric(1)))
  
  padding <- .1 * (fin - st) / n_bins
  bins <- seq(st - padding, fin + padding, length.out = n_bins + 1)
  return(bins)
}