
#' Plot Method for flirti Objects
#' 
#' Plots coefficient function of objects of class "flirti".
#'
#' @param x An object of class "flirti".
#' @param ... Other parameters to be passed through to plotting functions.
#'
#' @return A line graph of the beta values versus time.
#' @export
#'
#' @examples
#' library(fda)
#' betaind = 1
#' snr  = 2
#' nsim = 200
#' n    = 50
#' p    = 21
#' Y = array(NA,c(n,nsim))
#' X = array(NA,c(n,p,nsim))
#' domain = c(0,1)
#' lambda = seq(0.0005,0.01,length.out = 10)
#' Mf = 6:13
#' extra=list(Mf=Mf,lambda=lambda)
#' 
#' for(itersim in 1:nsim)
#' {
#'   dat = ngr.data.generator.bsplines(n=n,nknots=64,norder=4,p=p,domain=domain,snr=snr,betaind=1)
#'  Y[,itersim]  = dat$Y
#'   X[,,itersim] = dat$X
#' }
#' 
#' 
#' fltyfit = FLiRTI(Y=Y[1:n,1],(X[1:n,,1]),d=3,cons=4,domain=domain,extra=extra)
#' 
#' plot(fltyfit)
#' 
plot.flirti = function(x,...) {
  object  = x 
  plot(object$beta, type = "l", lwd = 2, ylim = c(min(object$beta), max(object$beta)), xlab = "Time", ylab = expression(beta(t)),...)
}

