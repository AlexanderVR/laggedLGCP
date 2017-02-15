#' Plot posterior quantiles of the latent rate functions.
#' 
#' @export
#' @param res Output list from lag_poisson_fit.
#' @param quant_level Which quantiles to use for shaded region.
#'
plot_results <- function(res, quant_level = 1) {
  
  quants <- res$quantiles$rates
  
  # index of 50% quantile
  med_idx <- ceiling(dim(quants)[1] / 2)
  
  # Select only the desired quantiles
  quants = quants[c(quant_level, med_idx, 2 * med_idx - quant_level), , ]
  par(mfrow=c(res$data$n_series, 1), mar=c(4,4,1,1))
  
  # Main title chould include the quantile names
  quantnames <- rownames(quants)
  maintitle = sprintf('Posterior median, %s and %s quantiles of posterior rates.', 
                      quantnames[1], quantnames[3])
  
  # append bin midpoints
  mids <- res$data$breaks[-1] - diff(res$data$breaks) / 2
  for (j in 1 : res$data$n_series) {
    # plot j'th process
    gp_plotting(mids, quants[,,j])
    if (res$data$coalescent) {
      if (j==1) {
        yl <- "Eff. Pop. Size"
        points(res$coal_data$coal_times, rep(0, length(res$coal_data$coal_times)), pch=3)
      } else {
        yl <- "Sampling Rate"
        points(res$coal_data$samp_times, rep(0, length(res$coal_data$samp_times)), pch=3)
      }
    } else {
      yl <- sprintf('Poisson rate for process %d', j)
      points(res$coal_data[[j]], rep(0, length(res$coal_data[[j]])), pch=3)
    }
    title(ylab=yl)
    if (j==1) title(main = maintitle)
    
  }
  title(xlab = res$data$time_unit)
}

gp_plotting <- function(midpoints, quantiles){
  
  temp_time = rep(midpoints, each=2)
  lower <- rep(quantiles[1, ], each=2)
  upper <- rep(quantiles[3, ], each=2)
  ymax <- max(quantiles[3, 1:(ceiling(dim(quantiles)[2] * .5))])
  x_points = c(temp_time, rev(temp_time))
  y_points = c(lower, rev(upper))

  # first plot medians
  plot(midpoints, quantiles[2, ], col="black",type="l",lwd=2.5, xlab="", ylab="", ylim=c(0, ymax))

  # plot interval defined by quantiles
  polygon(x_points, y_points, col=adjustcolor("gray",alpha.f=0.5),border="gray",lty=1,lwd=.5)
  points(midpoints,quantiles[1, ],col="black",type="l",lwd=0.5)
  points(midpoints,quantiles[3, ],col="black",type="l",lwd=0.5)
}

#' Plot posterior effective pop. size using phylodyn's plot_BNPR code
#' 
#' @export
#' 
#' @param fit Output of fit_LGCP_lag()
#' #param coal_data Original coalescent data -- a list containing samp_times, coal_times, and n_sampled elements.
#' @param traj True effective pop. size trajectory (if known)
#' @param ... Other args to pass to plot_BNPR
#' 
plot_coal_result <- function(fit, coal_data=NULL, traj = NULL, ...) {
  grid <- fit$data$breaks
  x <- grid[-1] - diff(grid) / 2
  effpop025 <- fit$quantiles$rates['2.5%', , 1]
  effpop975 <- fit$quantiles$rates['97.5%', , 1]
  effpop <- fit$quantiles$rates['50%', , 1]
  
  res <- c(coal_data, list(x=x, grid=grid, effpop=effpop, effpop025=effpop025, effpop975=effpop975))
  plot_BNPR(res, traj=traj, ...)
  
}
