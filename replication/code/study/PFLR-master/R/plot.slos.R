
#' Plot Method for SLoS Objects
#' 
#' Plots coefficient function for objects of class "slos".
#'
#' @import fda
#' @import psych
#' @import glmnet
#' @import stats
#' @import utils
#' @import MASS
#' @import flare
#' @param x An object of class "slos".
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
#' norder   = d+1
#' nknots   = M+1
#' knots    = seq(domain[1],domain[2],length.out = nknots)
#' nbasis   = nknots + norder - 2
#' basis    = create.bspline.basis(knots,nbasis,norder)
#' V = eval.penalty(basis,int2Lfd(2))
#'
#' extra=list(lambda=exp(seq(-18,-12, length.out = 10)),gamma=10^(-8:-6))
#'
#' for(itersim in 1:nsim)
#' {
#'  dat = ngr.data.generator.bsplines(n=n,nknots=64,norder=4,p=p,domain=domain,snr=snr,betaind=betaind)
#'   Y[,itersim]  = dat$Y
#'   X[,,itersim] = dat$X
#' }
#'
#' slosfit = SLoS(Y=Y[1:n,1],(X[1:n,,1]),M=M,d=d,domain=domain,extra=extra)
#' plot(slosfit)
#' 
#' 
plot.slos = function(x,...) {
  object = x  
  plot(object$beta, type = "l", lwd = 2, ylim = c(min(object$beta), max(object$beta)), xlab = "Time", ylab = expression(beta(t)),...)
}


