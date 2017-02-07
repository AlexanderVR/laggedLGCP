bin_coal <- function(coal_times, samp_times, n_sampled, n_bins) {
  # Approximation of Coalescent likelihood as a binned Poisson process
  #
  # Input: flu data of coalescent times, sampling times, number of lineages,
  # and number of grid points.
  #
  # Returns: list of counts y1, y2, mask for y2,

  # create grid of times
  all_times <- sort(c(samp_times, coal_times));
  grid <- seq(0, max(all_times), length.out = n_bins)
  midpts <- grid[-1] - diff(grid)/2

  # create coalescent likelihood for f_1 -- log effective pop size
  init <- coal_lik_init(samp_times=samp_times, n_sampled=n_sampled,
                        coal_times=coal_times, grid=grid)

  # Transform 'init' so that Coalescent(f_1) is Poisson with offset E
  # likelihood: Coalescent(f_1) = sum(y .* (-f_1) - E .* exp(-f_1))
  ind <- 1:(init$ng)

  # Create factors -- one per bin
  ind_rep = rep(ind, init$gridrep)

  # Coalescence over [t, t+delta] occurs at rate delta * l(l-1)/2, where
  # l = number of active lineages at time t.
  E_expanded = init$D * init$C

  # calculates the offset for each bin.
  E <- as.vector(tapply(E_expanded, ind_rep, sum))

  # There is a classic error that fails for zero elements of E,
  # so set those small values
  E <- ifelse(E==0, 1, E)

  # Number of coalescence events in each bin
  y <- as.vector(tapply(init$y, ind_rep, sum))

  # Calculate number of sampling events in each bin
  new_samp_times <- rep(0,length(grid)-1)
  for (j in 2:length(grid)){
    new_samp_times[(j-1)] <- sum(n_sampled[samp_times > grid[j-1] &
                                             samp_times <= grid[j]])
  }
  # add sampling time at boundary (t=0).
  new_samp_times[1] <- new_samp_times[1]+1

  # treat the counts (null) beyond the oldest sampling time as missing data
  miss <- sum(grid > max(samp_times))

  # mask for available sampling times data
  Nvis_samp <- init$ng - miss + 1

  # Grid size
  delta <- diff(grid)[1]

  # Compute the pairs of Poisson observations over grid points
  counts <- matrix(0, init$ng, 2)
  counts[, 1] <- y
  counts[1:Nvis_samp, 2] <- new_samp_times[1:Nvis_samp]

  offset <- matrix(0, init$ng, 2)
  offset[, 1] <- E
  offset[, 2] <- delta

  poisson_args <- list(counts=counts, breaks=grid, offset=offset, delta=delta,
                       N_vis=c(init$ng, Nvis_samp), N=init$ng, n_series=2, coalescent=1)

  return(poisson_args)

}
