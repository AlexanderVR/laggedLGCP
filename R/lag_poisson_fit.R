#' Fitting 1-D Log-Gaussian Cox Process with lag.
#' 
#' @export
#' 
#' @param events A list of vectors of timestamps, one vector for each observed process.
#' @param time_unit Unit of measurement for the timestamps.
#' @param likelihood Either 'poisson' or 'coalescent'. If 'coalescent', param EVENTS must contain
#' vectors named 'coal_times', 'samp_times', and 'n_sampled'.
#' @param n_bins Number of bins to use for binning the events
#' @param max_lag Place a Normal(0, sigma=max_lag/2) prior on the lag parameter.
#' @param max_n_frequencies Maximum number of cosine and sine basis functions to use in spectral representation.
#' @param min_percent_padding How far to extend the domain of the latent (periodic) Gaussian process.
#' @param prob_quantiles Which quantiles of the latent parameters should be returned?
#' @param return_stanfit If TRUE, the returned list includes the stanfit object containing all MCMC samples.
#' @param fitting_args Additional list of arguments that affect the marginal GMO fitting and MCMC simulation code.
#' 
#' @useDynLib laggedLGCP, .registration = TRUE
#'
fit_LGCP_lag <- function(events,
                            time_unit = 'Weeks',
                            likelihood = 'poisson',
                            n_bins = 100, 
                            max_lag = 5,
                            max_n_frequencies = 1023, 
                            min_percent_padding = 10.0, 
                            prob_quantiles = c(.025, .25),
                            return_stanfit = FALSE,
                            fitting_args = list()) {

  stopifnot(any(likelihood == c('poisson', 'coalescent')))
  
  # set defaults for GMO and stan fitting
  fitting_defaults <- list(n_chains = 1, 
                           n_samples = 1000,
                           refresh = 10, 
                           gmo_iter = 100, 
                           gmo_draws = 10, 
                           gmo_tol = 1e-4,
                           gmo_eta = 1.,
                           gmo_init_lengthscale = 10.,
                           gmo_max_block_size = 256)

  # Override defaults if specified
  for (arg_name in names(fitting_args)) {
    stopifnot(any(arg_name == names(fitting_defaults)))
    fitting_defaults[[arg_name]] <- fitting_args[[arg_name]]
  }
  
  stopifnot(max(prob_quantiles) < 0.5)
  
  # fill out probability quantiles with median and symmetric values
  prob_quantiles <- c(prob_quantiles, 0.5, rev(1 - prob_quantiles))
  
  # Make sure event data is a list of numeric vectors
  stopifnot(all(vapply(events, function(x) is.numeric(x) & is.vector(x), logical(1))))

  # Bin the event data
  if (likelihood == 'poisson') {
    pars <- bin_poisson(dat$events, n_bins=n_bins)
  } else {  # have coalescent data 
    pars <- bin_coal(events$coal_times, events$samp_times, events$n_sampled, n_bins)
  }
      
  # Set the prior on lag parameter to N(mu = 0, sigma = max_lag / 2)
  pars$max_lag <- max_lag

  # number of bins WITH padding -- must be a power of 2.
  pars$Ntot <- nextpow2(pars$N, min_percent_padding)

  # how many hyperparameters to optimize
  pars$npars <- 4  

  # number of frequencies to use
  pars$n_nonzero_freq <- min(max_n_frequencies, floor(pars$Ntot/2 - 1))

  # marginal process fitting function
  marginal_optimize <- function(k) {
    ind_pars <- split_data(pars, k)
    if (k==2 && pars$coalescent) {
      # second process (sampling) does not have the coalescent likelihood. 
      ind_pars$coalescent = 0
    }
    print(ind_pars$coalescent)
    sfit=sampling(stanmodels$single_gp_fixed_gmo, data = c(ind_pars, list(GMO_FLAG = FALSE, fixed_phi = double())),
         chains = 0, iter = 1)
    opt_res <- gmo(full_model=sfit, 
                 data = ind_pars, 
                 iter = as.integer(fitting_defaults$gmo_iter),
                 tol = fitting_defaults$gmo_tol, 
                 eta = fitting_defaults$gmo_eta, 
                 draws = as.integer(fitting_defaults$gmo_draws), 
                 init = get_init(ind_pars, lengthscale=fitting_defaults$gmo_init_lengthscale),
                 max_block_size = fitting_defaults$gmo_max_block_size)

    u_opt <- opt_res$par
    return(c(mu=u_opt[1], sigma=exp(u_opt[2]), lengthscale=exp(u_opt[3]), nu=exp(u_opt[4])))
  }
  # optimize hyperparameters of each marginal latent GP
  hyp_opt <- vapply(1:pars$n_series, FUN=marginal_optimize, FUN.VALUE=numeric(pars$npars))

  # Package the hyperparameters together
  pars <- c(pars, as.list(data.frame(t(hyp_opt))))

  # Now train joint GP stan model, with the above fixed hyperparameters
  fit <- sampling(stanmodels$multiple_gp_nonsep_fixed, 
                  data=pars, 
                  chains=fitting_defaults$n_chains, 
                  iter=fitting_defaults$n_samples, 
                  control=list(adapt_delta=0.95), 
                  refresh=fitting_defaults$refresh,
                  pars=c("re", "im", "x", "cor"), include=FALSE) # ignore high-dimensional latent variables

  # extract quantiles for all posterior samples
  samples <- rstan::extract(fit, permuted=TRUE)
  rates <- samples$rates
  quantiles <- list()
  for (param in names(samples)) {
    samp <- samples[[param]]
    if (is.vector(samp) || length(samp) == dim(samp)[1])
      quantiles[[param]] <- quantile(samp, probs = prob_quantiles)
    else {
      n_dim <- length(dim(samp))
      quantiles[[param]] <- apply(samp, MARGIN=2:n_dim,
                    FUN=function(x) quantile(x, probs=prob_quantiles))
    }
  }
  
  # Naive cross-correlation estimate of lag
  shifts_ccf <- rep(0, pars$n_series - 1)
  for (k in 1:(pars$n_series - 1)) {
    crosscor <- ccf(pars$counts[, pars$n_series], pars$counts[,k], plot=FALSE)
    shifts_ccf[k] <- crosscor$lag[which.max(crosscor$acf)] * pars$delta
  }

  # Data to return
  res <- list(quantiles = quantiles, shifts_ccf = shifts_ccf, data = pars, 
              events=events)
  if (return_stanfit)
    res$stanfit <- fit
  return(res)
}

#' Return marginal GP hyperparameters fitted by optimizing the marginal likelihood.
#' 
#' @export
#' @param res Result from fit_LGCP_lag() function.
#' 
get_hypers <- function(res)
  return(res$data[c('mu', 'sigma', 'lengthscale', 'nu')])
