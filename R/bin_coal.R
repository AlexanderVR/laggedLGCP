bin_coal <- function(coal_times, samp_times, n_sampled, n_bins, st=NULL, fin=NULL) {
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

  # extend grid if necessary
  ext = extend_grid(st, fin, grid)

  # Compute the pairs of Poisson observations over grid points
  counts <- matrix(-1, ext$n_grid, 2)
  counts[ext$old_idx, 1] <- y
  counts[ext$old_idx[1:Nvis_samp], 2] <- new_samp_times[1:Nvis_samp]

  offset <- matrix(-1, ext$n_grid, 2)
  offset[ext$old_idx, 1] <- E
  offset[ext$old_idx[1:Nvis_samp], 2] <- delta

  poisson_args <- list(counts=counts, breaks=ext$new_grid, offset=offset, delta=delta, N=ext$n_grid, n_series=2, coalescent=1)

  return(poisson_args)

}
#' Extend existing grid, keeping grid width constant.
#'
#' @param st Maximum start time of new grid
#' @param fin Minimum finish time of new grid
#' @param grid Existing grid to extend
extend_grid <- function(st, fin, grid) {
  dt = diff(grid)[1]
  if (is.null(st) || st == grid[1]) left <- c()
  else {
    stopifnot(st < 0)
    left = -rev(1 : ceiling(-st / dt)) * dt
  }

  if (is.null(fin) || fin == grid[length(grid)]) right <- c()
  else {
    stopifnot(fin > max(grid))
    right <- dt * (1 : ceiling(fin - max(grid)) / 2) + max(grid)
  }
  new_grid <- c(left, grid, right)
  idx <- (1 : (length(grid) - 1)) + length(left)
  return(list(left=left, right=right, new_grid=new_grid, old_idx=idx, n_grid=length(new_grid) - 1))
}