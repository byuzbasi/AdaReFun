# *** FLiRTI method (James et al 2009)

#############################
# *** FLiRTI Method
#############################
#---------------------------------------------#
# INPUT #
# y : vector of length n, centered response
# u : centered design matrix of n x p
# basis: bspline basis with degree d and M+1 equally spaced knots
# lambda: tuning parameter
# omega: tuning parameter
# d2: tuning parameter taking a value from 1,2,3,4
# nl: lenght.out of t observations that we need to calculate basisamt when estimating beta
# cons: divide subinterval into how many small ones
# flirty2 applies lasso as regularization method
#---------------------------------------------#
flirty1_dz2 = function(Y,U,basis,lambda,cons)
{

  n = length(Y) # number of observations
  p = dim(U)[2]

  # beta b-spline basis
  rng = getbasisrange(basis)
  breaks = c(rng[1],basis$params,rng[2])
  L = basis$nbasis
  M = length(breaks) - 1
  d = L-M

  t = seq(rng[1],rng[2],length.out=p)

  A = eval.basis(t, basis)

  V=U%*%solve(A)

  lenbeta = cons*p-3
  tobs=seq(rng[1],rng[2],length.out = lenbeta)
  basismat = eval.basis(tobs, basis)

  slimmm = slim(X=V, Y=Y, lambda = lambda, nlambda = 1,
                method="dantzig", q = 1, res.sd = FALSE,
                prec = 1e-7, max.ite = 1e5, verbose = FALSE)
  gamma = slimmm$beta
  eta = solve(A)%*%gamma
  beta=basismat%*%eta

  inv =t
  b=gamma
  b2=c(b[2:length(b)],0)
  b3 = abs(b)+abs(b2)
  inta = rep(NA,length(inv)-1)
  beta0 = rep(NA,lenbeta)
  for(aa in 1:length(inta))
  {
    if(b3[aa]==0){beta0[(cons*(aa-1)+1):(cons*aa+1)]=0}else{beta0[(cons*(aa-1)+1):(cons*aa+1)]=beta[(cons*(aa-1)+1):(cons*aa+1)]}
  }

  return(list(gamma = gamma,beta=beta0,A=A,df=slimmm$df))
}





#' FLiRTI Regression Model
#'
#' Calculates functional regression that's interpretable using the FLiRTI method.
#'
#' @import fda
#' @import psych
#' @import glmnet
#' @import stats
#' @import utils
#' @import MASS
#' @import flare
#' @param Y Vector of length n, centred response.
#' @param X Matrix of n x p, covariate matrix, should be dense.
#' @param d Integer, degree of the B-spline basis functions.
#' @param cons Divide subinterval into how many small ones.
#' @param domain The range over which the function X(t) is evaluated and the coefficient function \eqn{\beta}(t) is expanded by the B-spline basis functions.
#' @param extra List containing parameters which have default values:
#'              \itemize{
#'                \item Mf: Mf+1 is the number of knots for the B-spline basis functions that expand \eqn{\beta}(t), default is 6:30.
#'                \item lambda: Tuning parameter, default is seq(0.0005,100,length.out = 50).
#'              }
#'
#' @return beta: Estimated \eqn{\beta}(t) at discrete points.
#' @return extra: List containing other values which may be of use:
#'        \itemize{
#'          \item X: Matrix of n x p used for model.
#'          \item Y: Vector of length n used for model.
#'          \item domain: The range over which the function X(t) was evaluated and the coefficient function \eqn{\beta}(t) was expanded by the B-spline basis functions.
#'          \item delta: Estimated cutoff point.
#'          \item OptM: Optimal number of B-spline knots selected by BIC.
#'          \item Optlambda: Optimal shrinkage parameter selected by BIC.
#'        }
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
#'
FLiRTI = function(Y,X,d,cons,domain,extra=list(Mf=6:30,lambda=seq(0.0005,100,length.out=50)))
{
  X = t(X)
  # lambda is a sequence of values
  df = d
  if (is.null(extra$Mf)) {
    Mf=6:30
  } else {
    Mf=extra$Mf
  }

  if (is.null(extra$lambda)) {
    lambda=seq(0.0005,100,length.out=50)
  } else {
    lambda=extra$lambda
  }

  n = length(Y)
  p=dim(X)[1]
  tobs = seq(domain[1],domain[2],length.out = p)

  h = (domain[2]-domain[1])/(p-1)
  cef = c(1, rep(c(4,2), (p-3)/2), 4, 1)

  beta_flty = array(NA,c(100,length(Mf),length(lambda)))
  BICdant = array(NA,c(length(Mf),length(lambda)))


  for(i in 1:length(Mf))
  {
    nknots   = Mf[i] + 1
    knots    = seq(domain[1],domain[2],length.out=nknots)
    norder   = df + 1
    nbasis   = nknots + norder - 2
    basis    = create.bspline.basis(knots,nbasis,norder)
    basismat = eval.basis(tobs, basis)
    nbetat   = cons*(nbasis-1)+1

    U = h/3*t(X)%*%diag(cef)%*%basismat

    for(j in 1:length(lambda))
    {
      flirtyfit = flirty1_dz2(Y=Y,U=U,basis=basis,lambda=lambda[j],cons=cons)
      beta_flty[1:nbetat,i,j] = flirtyfit$beta
      A = flirtyfit$A

      dfdantzig = flirtyfit$df
      Yhat=U%*%solve(A)%*%flirtyfit$gamma
      res = sum((Y-Yhat)^2)/n
      BICdant[i,j] = n*log(res)+dfdantzig*log(n)
    }
  }
  idx = which(BICdant==min(BICdant,na.rm = T),arr.ind = T)
  OptM = Mf[idx[1]]
  Optlambda = lambda[idx[2]]

  Optnbasis   = OptM + norder - 1
  nbetat = cons*(Optnbasis-1)+1
  betaBIC = beta_flty[1:nbetat,idx[1],idx[2]]
  T1 = T1.hat(betaBIC)

  result = list(beta=betaBIC,extra=list(X=X, Y = Y, domain = domain, delta=T1,OptM=OptM,Optlambda=Optlambda))
  class(result) = 'flirti'
  result
}
