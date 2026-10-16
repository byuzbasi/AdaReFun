
#' Predict Method for Penalized B-splines objects
#'
#' Predicted values based on objects of class "ps".
#'
#' @import fda
#' @import psych
#' @import glmnet
#' @import stats
#' @import utils
#' @import MASS
#' @import flare
#' @param object An object of class "ps".
#' @param Xnew New covariate matrix for prediction, should be dense, centred.
#' @param ... Not applicable
#'
#' @return Predicted values.
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
#' predict(psfit,X[1:n,,1])
#'
predict.ps = function(object,Xnew,...) {
  psobj = object
  # Xnew should be a matrix
  nxnew = dim(Xnew)[1]
  pxnew = dim(Xnew)[2]
  Xmean = psobj$extra$Xmean
  Ymean = psobj$extra$Ymean
  xc = Xnew - matrix(rep(Xmean,nxnew),nrow=nxnew,byrow=TRUE)
  Ucnew = compute.u_PS(X=xc, M=psobj$extra$M, d=psobj$extra$d, domain=psobj$extra$domain)$U
  bps = psobj$extra$b
  Yhat = Ymean + Ucnew%*%bps
  return(Yhat)
}
