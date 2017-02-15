#' Fit a binned 1-D Log-Gaussian Cox Process with lag.
#' 
#' @param binned_data A list of type BinnedData. Usually obtained from bin_coalescent_data() or bin_poisson_data().
#' @param max_lag Lag parameter constrained to lie within [-max_lag, max_lag].
#' @param prior_lengthscale_mean Where to center the LogNormal(., 1) prior on the lengthscale parameter
#' @param prior_smoothness_mean Mean of Normal prior for the smoothness parameter 'nu'.
#' @param prior_smoothness_std Standard deviation of Normal prior for the smoothness parameter 'nu'.
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
#' # Bin the events then fit the LGCP
#' binned_data <- bin_coalescent_data(sim, n_bins=50)
#' fit <- fit_binned_LGCP(binned_data, max_lag=3)
#' # Plot inferred and 'true' effective population size trajectory
#' par(mfrow=c(3,1))
#' plot_coal_result(fit, coal_data= sim$coal_data, traj=sim$rate_functions$coal_fun, main="Lagged Preferential Sampling", ylim=c(0.5,15))
#' # Plot non-Pref. sampling and PS without lag inferences``
#' plot_BNPR(BNPR(sim$coal_data, lengthout = 100), traj=sim$rate_functions$coal_fun, main="No PS", ylim=c(0.5,15))
#' plot_BNPR(BNPR_PS(sim$coal_data, lengthout = 100), traj=sim$rate_functions$coal_fun, main="PS without lag", ylim=c(0.5,15))
#' @export
fit_binned_LGCP <- function(binned_data,
                            max_lag,
                            prior_lengthscale_mean = 1.0,
                            prior_smoothness_mean = 1.0,
                            prior_smoothness_std = 0.5,
                            max_n_frequencies = 1023, 
                            min_percent_padding = 10.0, 
                            prob_quantiles = c(.025, .25),
                            return_stanfit = FALSE,
                            fitting_args = list()) {

  stopifnot(class(binned_data) == "BinnedData")
  stopifnot(max_lag > 0)
  stopifnot(max(prob_quantiles) < 0.5)
  
  # set defaults for GMO and stan fitting
  fitting_defaults <- list(n_chains = 1, 
                           n_samples = 1000,
                           refresh = 10, 
                           gmo_iter = 100, 
                           gmo_draws = 10, 
                           gmo_tol = 1e-4,
                           gmo_eta = .5,
                           gmo_max_block_size = 256)

  # Override defaults if specified
  for (arg_name in names(fitting_args)) {
    stopifnot(any(arg_name == names(fitting_defaults)))
    fitting_defaults[[arg_name]] <- fitting_args[[arg_name]]
  }
  
  # add priors to data to send to stan
  binned_data = c(binned_data, list(prior_lengthscale_mean=prior_lengthscale_mean,
                                    prior_smoothness_mean=prior_smoothness_mean,
                                    prior_smoothness_std=prior_smoothness_std))

  # fill out probability quantiles with median and symmetric values
  prob_quantiles <- c(prob_quantiles, 0.5, rev(1 - prob_quantiles))

  # Set the prior on lag parameter to N(mu = 0, sigma = max_lag / 2)
  binned_data$max_lag <- max_lag

  # number of bins WITH padding -- must be a power of 2.
  binned_data$Ntot <- nextpow2(binned_data$N, min_percent_padding)

  # how many hyperparameters to optimize
  binned_data$npars <- 4  

  # number of frequencies to use
  binned_data$n_nonzero_freq <- min(max_n_frequencies, floor(binned_data$Ntot/2 - 1))

  # marginal process fitting function
  marginal_optimize <- function(k) {
    ind_pars <- split_data(binned_data, k)
    
    if (k > 1 && binned_data$coalescent) {
      # second process (sampling) does not have the coalescent likelihood. 
      ind_pars$coalescent = 0
    }
    
    sfit = sampling(stanmodels$single_gp_fixed_gmo, data = c(ind_pars, list(GMO_FLAG = FALSE, fixed_phi = double())), chains = 0, iter = 1)

    opt_res <- gmo(full_model=sfit, 
                 data = ind_pars, 
                 iter = as.integer(fitting_defaults$gmo_iter),
                 tol = fitting_defaults$gmo_tol, 
                 eta = fitting_defaults$gmo_eta, 
                 draws = as.integer(fitting_defaults$gmo_draws), 
                 init = get_init(ind_pars, 
                                 lengthscale=prior_lengthscale_mean,
                                 nu=prior_smoothness_mean),
                 max_block_size = fitting_defaults$gmo_max_block_size)

    u_opt <- opt_res$par
    return(c(mu=u_opt[1], sigma=exp(u_opt[2]), lengthscale=exp(u_opt[3]), nu=exp(u_opt[4])))
  }
  # optimize hyperparameters of each marginal latent GP
  hyp_opt <- vapply(1:binned_data$n_series, FUN=marginal_optimize, FUN.VALUE=numeric(binned_data$npars))

  # Package the hyperparameters together
  binned_data <- c(binned_data, list(N_vis = sum(binned_data$counts >= 0)), as.list(data.frame(t(hyp_opt))))

  # Now train joint GP stan model, with the above fixed hyperparameters
  fit <- sampling(stanmodels$multiple_gp_nonsep_fixed, 
                  data=binned_data, 
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
  shifts_ccf <- rep(0, binned_data$n_series - 1)
  for (k in 1:(binned_data$n_series - 1)) {
    crosscor <- ccf(binned_data$counts[, binned_data$n_series], binned_data$counts[,k], plot=FALSE)
    shifts_ccf[k] <- crosscor$lag[which.max(crosscor$acf)] * binned_data$delta
  }

  # Data to return
  res <- list(quantiles = quantiles, shifts_ccf = shifts_ccf, data = binned_data)
  if (return_stanfit)
    res$stanfit <- fit
  return(res)
}

#' Return marginal GP hyperparameters fitted by optimizing the marginal likelihood.
#' 
#' @param res Result from fit_binned_LGCP() function.
#'
#' @export
get_hypers <- function(res)
  return(res$data[c('mu', 'sigma', 'lengthscale', 'nu')])
