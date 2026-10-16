#' Summary Method for Penalized B-splines Objects
#'
#' Summarizes the values of an object of class "ps".
#'
#' @param object An object of class "ps".
#' @param ... Not applicable
#'
#' @return Prints a 5 number summary of the beta values and coefficient values, and the optimal alpha.
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
#' summary(psfit)
#'
summary.ps = function(object,...) {
  x=object
  print("Summary of estimated beta values:")
  cat("\n")
  print(summary(x$beta))
  cat("\n")
  print("Summary of estimated B-spline coefficient values:")
  cat("\n")
  print(summary(x$extra$b))
  cat("\n")
  print(paste0("Optimal alpha value chosen: ",x$extra$Optalpha))
}
