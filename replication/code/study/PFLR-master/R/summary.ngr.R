

#' Summary Method for ngr Objects
#' 
#' Summarizes objects of class "ngr".
#'
#' @param object An object of class "ngr".
#' @param ... Not applicable
#'
#' @return Prints the 5 number summaries of beta and b values. Prints delta, Optkappa, and Opttau values.
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
  #' dat = ngr.data.generator.bsplines(n=n,nknots=64,norder=4,p=p,domain=domain,snr=snr,betaind=1)
  #' Y[,itersim]  = dat$Y
  #' X[,,itersim] = dat$X
#' }
#'
#' ngrfit = ngr(Y=Y[1:n,1],X=(X[1:n,,1]),M,d,domain,extra= list(alphaPS=alphaPS, kappa=kappa, tau=tau))
#' 
#' summary(ngrfit)
#' 
summary.ngr = function(object,...)
{
  ngrobj = object
  print("Summary of estimated beta values:")
  cat("\n")
  print(summary(ngrobj$beta))
  cat("\n")
  print("Summary of estimated b values:")
  cat("\n")
  print(summary(ngrobj$extra$b))
  cat("\n")
  print(paste0("Estimated cutoff point: ",ngrobj$extra$delta))
  cat("\n")
  print(paste0("Optimal roughness penalty selected: ",ngrobj$extra$Optkappa))
  cat("\n")
  print(paste0("Optimal group bridge penalty selected: ",ngrobj$extra$Opttau))
}

