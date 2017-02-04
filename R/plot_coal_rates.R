#' Plot posterior quantiles of the latent rate functions.
#' 
#' @export
#' @param res Output list from lag_poisson_fit.
#'
plot_coal_results <- function(res) {

  # extract posterior rates
  pars <- res$pars
  samples <- rstan::extract(res$fit, permuted=TRUE)
  rates <- samples$rates
  d <- dim(rates)[3]
  par(mfrow=c(d,1), mar=c(4,4,1,1))
  # append bin midpoints
  mids <- (1:dim(rates)[2]) * pars$delta
  for (j in 1:d) {
    rj <- rates[,,j]
    if (j==1)
      rj = 1/rj
    # get point-wise quantiles of posterior rates
    quants <- plyr::aaply(rj, .margins=2,
                    .fun=function(x) quantile(x, probs=c(0.05, 0.5, 0.95)))
    # plot results
    gp_plotting(mids, quants)
    if (j==1)
      yl <- "Eff. Pop. Size"
    else
      yl <- "Sampling Rate"
    title(ylab=yl, xlab = "Week")
    #points(coal$date, rep(1, length(coal$date)), pch=3)
  }
}

gp_plotting <- function(midpoints, quantiles){
  temp_time = rep(midpoints, each=2)
  lower <- rep(quantiles[,1], each=2)
  upper <- rep(quantiles[,3], each=2)
  ymax <- max(quantiles[1:(ceiling(dim(quantiles)[1] * .5)),3])
  x_points = c(temp_time, rev(temp_time))
  y_points = c(lower, rev(upper))

  # first plot medians
  plot(midpoints, quantiles[,2], col="black",type="l",lwd=2.5, xlab="", ylab="", ylim=c(0, ymax))

  # plot interval defined by quantiles
  polygon(x_points, y_points, col=adjustcolor("gray",alpha.f=0.5),border="gray",lty=1,lwd=.5)
  points(midpoints,quantiles[,1],col="black",type="l",lwd=0.5)
  points(midpoints,quantiles[,3],col="black",type="l",lwd=0.5)
}
