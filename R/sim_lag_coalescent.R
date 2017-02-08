#' Generate a realization of a lagged Coalescent process. 
#' 
#' Returns a list of sampling times, coalescent times, and n_samples.
#' 
#' @export
#' 
#' @param lag The 'true' lag parameter, indicating how much the effective population size
#'  N_e(t) is shifted to the right relative to the sampling rate f(t).
#' relative to the last series.
#' @param type Which type of latent trajectory to use. Either 'cyclic' or 'boombust'.
#' @param c Constant parameter in equation \eqn{f(t) = c N_e(t)^\beta}.
#' @param beta Power parameter in above equation relating sampling rate f(t) to 
#' effective population size N_e(t).
#' @param scaling Constant scaling factor for the overall Poisson rate.
#' 
#' @importFrom phylodyn cyclic_traj boombust_traj pref_sample coalsim BNPR BNPR_PS
#' 
sim_lag_coalescent <- function(lag=0, type='cyclic', c=50, beta=0.5, scaling=.1) {
  
  if (type == 'cyclic')
    coal_fun <- function(t) cyclic_traj(t) * scaling
  else if (type == 'boombust')
    coal_fun <- function(t) boombust_traj(t, bust=15) * scaling
  else stop("unrecognized trajectory type")
  
  # lagged sampling rate function
  samp_fun <- function(x) coal_fun(x - lag)
  
  # generate sampling times from fun
  sampling_times <- pref_sample(samp_fun, lim=c(0, 40), c=c, beta=beta, upper=NULL, grid.len=2000)
  n_samples <- c(3, rep(1, length(sampling_times)))
  sampling_times <- c(0, sampling_times)
  
  # now generate coalescent times given sampling times and lagged N_e(t) trajectory
  sim <- coalsim(samp_times=sampling_times, n_sampled=n_samples, traj=coal_fun)
  events = sim[c('coal_times', 'samp_times', 'n_sampled')]
  
  return(list(events=events, lag=lag, type=type,
              rate_functions=list(coal_fun=coal_fun, samp_fun=samp_fun)))
}