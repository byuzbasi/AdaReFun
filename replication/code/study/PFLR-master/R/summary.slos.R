
#' Summary Method for SLoS Objects
#'
#' Summarizes values of an object of class "slos".
#'
#' @import fda
#' @import psych
#' @import glmnet
#' @import stats
#' @import utils
#' @import MASS
#' @import flare
#' @param object An object of class "slos".
#' @param ... Not applicable
#'
#' @return Prints five number summary of beta values, delta, Optgamma, and Optlambda.
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
#'
#' summary(slosfit)
#'
#'
summary.slos = function(object,...) {
  x = object 
  print("Summary of estimated beta values:")
  cat("\n")
  print(summary(x$beta))
  cat("\n")
  print(paste0("Estimated cutoff point: ",x$extra$delta))
  cat("\n")
  print(paste0("Optimal smoothing parameter: ",x$extra$Optgamma))
  cat("\n")
  print(paste0("Optimal shrinkage parameter: ",x$extra$Optlambda))
}

