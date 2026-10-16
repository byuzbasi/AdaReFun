
#' Predict Method for flirti Objects
#'
#' Predicted values based on objects of the class "flirti".
#'
#' @param object An object of class "flirti".
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
#' predict(fltyfit,(X[1:n,,1]))
#'
#'
predict.flirti = function(object,Xnew,...) {
  flobj = object
  domain = flobj$extra$domain
  # transfer beta into an fdobject
  tobs = seq(domain[1], domain[2], length.out = length(flobj$beta))
  smoothbeta = smooth.spline(tobs, flobj$beta, spar = 10^(-8))
  tt = seq(domain[1], domain[2], length.out = dim(flobj$extra$X)[1])
  # Smooth for covariates, centered
  betadense = predict(smoothbeta, tt)

  nxnew = dim(Xnew)[1]
  pxnew = dim(Xnew)[2]
  Xmean = apply(flobj$extra$X, 1, mean)
  Ymean = mean(flobj$extra$Y)
  xc = Xnew - matrix(rep(Xmean,nxnew),nrow=nxnew,byrow=TRUE)

  h   = (domain[2]-domain[1])/(length(tt)-1)
  Yhat= Ymean + h*xc%*%as.matrix(betadense$y)
  return(Yhat)
}
