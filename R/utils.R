
nextpow2 <- function(n, perc=10) {
  n <- n * (1+perc/100)
  return(2^ceiling(log2(n)))
}

split_data <- function(pars, j) {
  # return pars for the j'th univariate series, given pars for all.

  new_pars <- pars
  new_pars$offset <- pars$offset[, j]
  new_pars$N_vis <- pars$N_vis[j]
  new_pars$counts <- pars$counts[, j]
  return(new_pars)
}

get_pars <- function(fit) {
  opts <- rstan::optimizing(fit$local_model, data=c(fit$data, list(fixed_phi=fit$par)),
                          as_vector=FALSE)
  return(opts$par)
}

get_joint_init <- function(pars1, pars2, lengthscale=1, nu=1) {
  mu <- c(get_init(pars1)[1], get_init(pars2)[1])
  sigma <- c(get_init(pars1)[2], get_init(pars2)[2])
  corr = 0.
  shifts=0.
  return(c(mu, sigma, shifts, corr, log(lengthscale), log(nu)))
}

get_init <- function(pars, lengthscale = 1, nu = 1) {
  # return good initial values for GP hyperparameters
  mu <- log(0.00001 + mean(pars$counts[1:pars$N_vis] / pars$offset[1:pars$N_vis]))
  sigma <- 1.
  return(c(mu, log(sigma), log(lengthscale), log(nu)))
}

find_contiguous <- function(x, minbreak, st=0, fin=Inf) {
  # find indices of x with adjacent values at least 'minbreak' apart
  # return start and end values of contiguous segments
  ix <- which(x[-1] - x[1:(length(x) - 1)] > minbreak)
  if (fin == Inf) fin <- max(x)
  breaks <- sort(c(st, x[ix], x[ix + 1], fin))
  breaks <- matrix(breaks, 2)
  rownames(breaks) <- c("start", "end")
  return(breaks)
}

