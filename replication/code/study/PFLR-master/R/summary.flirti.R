

#' Summary Method for flirti Objects
#' 
#' Summarizes the values of an object of class "flirti".
#'
#' @param object An object of class "flirti".
#' @param ... Not applicable
#'
#' @return Prints a 5 number summary of the beta values, delta, OptM, and Optlambda
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
#' summary(fltyfit)
#' 
#' 
summary.flirti = function(object,...) {
  x = object 
  print("Summary of estimated beta values:")
  cat("\n")
  print(summary(x$beta))
  cat("\n")
  print(paste0("Estimated cutoff point: ",x$extra$delta))
  cat("\n")
  print(paste0("Optimal number of B-spline knots: ",x$extra$OptM))
  cat("\n")
  print(paste0("Optimal shrinkage parameter: ",x$extra$Optlambda))
}


