
#' Predict Method for SLoS objects
#'
#' Predicted values based on objects of class "slos".
#'
#' @import fda
#' @import psych
#' @import glmnet
#' @import stats
#' @import utils
#' @import MASS
#' @import flare
#' @param object An object of class "slos".
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
#' predict(slosfit,(X[1:n,,1]))
#'
predict.slos = function(object,Xnew,...) {
  slosobj = object 
  # Xnew should be a matrix
  nxnew = dim(Xnew)[1]
  pxnew = dim(Xnew)[2]
  Xmean = apply(slosobj$extra$X, 1, mean)
  Ymean = mean(slosobj$extra$Y)
  xc = Xnew - matrix(rep(Xmean,nxnew),nrow=nxnew,byrow=TRUE)
  Ucnew = compute.u(X=xc, n=nxnew, M=slosobj$extra$M, d=slosobj$extra$d, domain=slosobj$extra$domain)$U
  bslos = slosobj$extra$b
  Yhat = Ymean + Ucnew%*%bslos
  return(Yhat)
}

