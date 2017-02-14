#' Fitting 1-D Log-Gaussian Coalescent Process with lag and optional covariates
#' 
#' @export
#' 
#' @param coal_data A list containing elements 'coal_times', 'samp_times', and 'n_sampled'
#' @param max_lag Lag parameter constrained to lie within [-max_lag, max_lag].
#' @param covariate If specified, model lag between eff. pop. size and covariate, instead of sampling rate.
#' @param time_unit Unit of measurement for the timestamps.
#' @param n_bins Number of bins to use for binning the events.
#' @param max_n_frequencies Maximum number of cosine and sine basis functions to use in spectral representation.
#' @param min_percent_padding How far to extend the domain of the latent (periodic) Gaussian process.
#' @param prob_quantiles Which quantiles of the latent parameters should be returned?
#' @param return_stanfit If TRUE, the returned list includes the stanfit object containing all MCMC samples.
#' @param fitting_args Additional list of arguments that affect the marginal GMO fitting and MCMC simulation code.
#' 
#' @useDynLib laggedLGCP, .registration = TRUE
#'
#' @examples 
#' # To simulate a coalescent process with lagged sampling times
#' sim <- sim_lag_coalescent(lag=-1, c=1, beta=2, scaling=0.05)
#' fit <- fit_LGCP_lag(sim$events, max_lag=3, likelihood='coalescent')
#' # Plot inferred and 'true' effective population size trajectory
#' par(mfrow=c(3,1))
#' plot_coal_result(fit, traj=sim$rate_functions$coal_fun, main="Lagged Preferential Sampling", ylim=c(0.5,15))
#' # Plot non-Pref. sampling and PS without lag inferences``
#' plot_BNPR(BNPR(sim$events, lengthout = 100), traj=sim$rate_functions$coal_fun, main="No PS", ylim=c(0.5,15))
#' plot_BNPR(BNPR_PS(sim$events, lengthout = 100), traj=sim$rate_functions$coal_fun, main="PS without lag", ylim=c(0.5,15))

lag_coal_fit <- function(events,
                         max_lag,
                         covariate = NULL,
                         time_unit = 'Weeks',
                         n_bins = 100, 
                         max_n_frequencies = 1023, 
                         min_percent_padding = 10.0, 
                         prob_quantiles = c(.025, .25),
                         return_stanfit = FALSE,
                         fitting_args = list()) {

  # if (all(sort(names(events)) == c('coal_times', 'n_sampled', 'samp_times')))
  #   likelihood <- 'coalescent'
  # stopifnot(any(likelihood == c('poisson', 'coalescent')))
  stopifnot(max_lag > 0)
  stopifnot(max(prob_quantiles) < 0.5)
  
  # set defaults for GMO and stan fitting
  fitting_defaults <- list(n_chains = 1, 
                           n_samples = 1000,
                           refresh = 10, 
                           gmo_iter = 100, 
                           gmo_draws = 15, 
                           gmo_tol = 1e-4,
                           gmo_eta = 1.,
                           gmo_init_lengthscale = -1,
                           gmo_max_block_size = 256)

  # Override defaults if specified
  for (arg_name in names(fitting_args)) {
    stopifnot(any(arg_name == names(fitting_defaults)))
    fitting_defaults[[arg_name]] <- fitting_args[[arg_name]]
  }
  
  # adjust starting lengthscale value if unspecified
  if (fitting_defaults$gmo_init_lengthscale <= 0) {
    fitting_defaults$gmo_init_lengthscale <- 2 * max_lag
  }
  
  # fill out probability quantiles with median and symmetric values
  prob_quantiles <- c(prob_quantiles, 0.5, rev(1 - prob_quantiles))
  
  # Make sure event data is a list of numeric vectors
  stopifnot(all(vapply(events, function(x) is.numeric(x) & is.vector(x), logical(1))))

  # if have covariate extend the grid to accomodate
  if (!is.null(covariate)) {
    st <- min(c(covariate, 0))
    fin <- max(c(covariate, events$samp_times, events$coal_times))
  } else {
    st <- NULL
    fin <- NULL
  }

  # Bin the event data
  pars <- bin_coal(events$coal_times, events$samp_times, events$n_sampled, n_bins, st=st, fin=fin)

  # replace sampling times with covariate if specified
  if (!is.null(covariate)) {
    pars_cov <- bin_poisson(covariate, bins=pars$breaks)
    pars$counts <- cbind(pars$counts[, 1], pars_cov$counts)
    pars$offset <- cbind(pars$offset[, 1], pars_cov$offset) 
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
  pars <- c(pars, list(N_vis = sum(pars$counts >= 0)), as.list(data.frame(t(hyp_opt))))

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
  print(quantiles$shifts)
  
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
