
#' Predict Method for ngr Objects
#'
#' Predicted values based on "ngr" class objects.
#'
#' @param object An object of class "ngr".
#' @param Xnew New covariate matrix for prediction, should be dense, centred.
#' @param ... Not applicable
#'
#' @return Estimated Y hat value.
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
#' tobs = seq(domain[1],domain[2],length.out = p)
#' knots    = seq(domain[1],domain[2],length.out = nknots)
#' nbasis   = nknots + norder - 2
#' basis    = create.bspline.basis(knots,nbasis,norder)
#' basismat = eval.basis(tobs, basis) 
#' h = (domain[2]-domain[1])/M
#' cef = c(1, rep(c(4,2), (M-2)/2), 4, 1)
#'
#' V = eval.penalty(basis,int2Lfd(2))
#' alphaPS = 10^(-(10:3))
#' kappa   = 10^(-(8:7))
#' tau     = exp(seq(-35,-28,len=20))
#' gamma   = 0.5
#'
#'
#' for(itersim in 1:nsim)
#' {
  #' dat = ngr.data.generator.bsplines(n=n,nknots=64,norder=4,p=p,domain=domain,snr=snr,betaind=betaind)
  #' Y[,itersim]  = dat$Y
  #' X[,,itersim] = dat$X
#' }
#'
#' ngrfit = ngr(Y=Y[1:n,1],X=(X[1:n,,1]),M,d,domain,extra= list(alphaPS=alphaPS, kappa=kappa, tau=tau))
#' predict(ngrfit,X[1:n,,1])
#'
predict.ngr = function(object, Xnew,...)
{
  ngrobj = object 
  # Xnew should be a matrix
  nxnew = dim(Xnew)[1]
  pxnew = dim(Xnew)[2]
  xc = Xnew - matrix(rep(ngrobj$extra$Xmean,nxnew),nrow=nxnew,byrow=TRUE)
  Ucnew = compute.u(X=xc, n=nxnew, M=ngrobj$extra$M, d=ngrobj$extra$d, domain=ngrobj$extra$domain)$U
  bngr = ngrobj$extra$b
  Yhat = ngrobj$extra$Ymean + Ucnew%*%bngr
  return(Yhat)
}
