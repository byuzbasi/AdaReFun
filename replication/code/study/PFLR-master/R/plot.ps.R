#' Plot Method for Penalized B-splines Objects
#'
#' Plots coefficient function of objects of class "ps".
#'
#' @param x An object of class "ps".
#' @param ... Other parameters to be passed through to plotting functions.
#'
#' @return A line graph of the beta values versus time.
#' @export
#'
#' @examples
#' library(fda)
#' betaind = 1
#' snr  = 2
#' nsim = 1
#' n    = 50
#' p    = 21
#' Y = array(NA,c(n,nsim))
#' X = array(NA,c(n,p,nsim))
#' domain = c(0,1)
#'
#' M = 20
#' d = 3
#' alpha = 10^(-(10:3))
#'
#'
#' for(itersim in 1:nsim)
#' {
#' dat = ngr.data.generator.bsplines(n=n,nknots=64,norder=4,p=p,domain=domain,snr=snr,betaind=betaind)
#' Y[,itersim]  = dat$Y
#' X[,,itersim] = dat$X
#' }
#'
#' psfit = PenS(Y=Y[1:n,1],X=(X[1:n,,1]), alpha=alpha, M=M, d=d, domain=domain)
#'
#' plot(psfit)
#'
plot.ps = function(x,...) {
  object  = x
  plot(object$beta, type = "l", lwd = 2, ylim = c(min(object$beta), max(object$beta)), xlab = "Time", ylab = expression(beta(t)),...)
}
