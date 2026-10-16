#library(psych)
#' SLoS regression Model
#'
#' Calculates functional regression using the Smooth and Locally Sparse (SLoS) method.
#'
#' @import fda
#' @import psych
#' @import glmnet
#' @import stats
#' @import utils
#' @import MASS
#' @import flare
#' @param Y Vector, length n, centred response.
#' @param X Matrix of n x p, covariate matrix, should be dense, centred.
#' @param M Integer, t1,..., tM are M equally spaced knots.
#' @param d Integer, the degree of B-Splines.
#' @param domain The range over which the function X(t) is evaluated and the coefficient function \eqn{\beta}(t) is expanded by the B-spline basis functions.
#' @param extra List of parameters which have default values:
#'              \itemize{
#'                \item Maxiter: Maximum number of iterations for convergence of beta, default is 100.
#'                \item lambda: Positive number, tuning parameter for fSCAD penalty, default is exp(seq(-20,-12, length.out = 10)).
#'                \item gamma: Positive number, tuning parameter for the roughness penalty, default is 10^(-9:0).
#'                \item absTol: Number, if max(norm(bHat)) is smaller than absTol, we stop another iteration, default is 10^(-10).
#'                \item Cutoff: Number, if bHat is smaller than Cutoff, set it to zero to avoid being numerically unstable, default is 10^(-6).
#'              }
#'
#' @return beta: Estimated \eqn{\beta}(t) at discrete points.
#' @return extra: List containing other values which may be of use:
#'        \itemize{
#'          \item X: Matrix of n x p used for model.
#'          \item Y: Vector of length n used for model.
#'          \item M: Integer representing the number of knots used in the model calculation.
#'          \item d: Integer, degree of B-Splines used.
#'          \item domain: The range over which the function X(t) was evaluated and the coefficient function \eqn{\beta}(t) was expanded by the B-spline basis functions.
#'          \item b: Estimated b values.
#'          \item delta: Estimated cutoff point.
#'          \item Optgamma: Optimal smoothing parameter selected by BIC.
#'          \item Optlambda: Optimal shrinkage parameter selected by BIC.}
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
#' basis    = fda::create.bspline.basis(knots,nbasis,norder)
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
#'
SLoS = function(Y,X,M,d,domain,extra=list(Maxiter=100,lambda=exp(seq(-20,-12, length.out = 10)),gamma=10^(-9:0),absTol=10^(-10),Cutoff=10^(-6)))
{
 X = t(X)
  if (is.null(extra$Maxiter)) {
    Maxiter = 100
  } else {
    Maxiter=extra$Maxiter
  }

  if (is.null(extra$lambda)) {
    lambda=exp(seq(-30,0,length.out=10))
  } else {
    lambda=extra$lambda
  }

  if (is.null(extra$gamma)) {
    gamma=10^(-10:10)
  } else {
    gamma=extra$gamma
  }

  if (is.null(extra$absTol)) {
    absTol=10^(-10)
  } else {
    absTol=extra$absTol
  }

  if (is.null(extra$Cutoff)) {
    Cutoff=10^(-6)
  } else {
    Cutoff=extra$Cutoff
  }

  Yn = Y
  Xn = t(X)
  Ymean=mean(Yn)
  Xmean = apply(Xn,2,mean) # column mean of the covariate curves
  n = length(Y)
  m = dim(X)[2]

  norder   = d+1
  knots    = seq(domain[1],domain[2], length.out=M+1)
  nknots   = length(knots)
  nbasis   = length(knots) + norder - 2 # i+2
  beta.basis  = create.bspline.basis(knots,nbasis,norder)
  V = eval.penalty(beta.basis,int2Lfd(2))
  beta.slos.all = array(NA,c(M+1,length(lambda),length(gamma)))
  b.slos.all = array(NA,c(M+d,length(lambda),length(gamma)))
  BIC_slos  = array(NA,c(length(lambda),length(gamma)))

  #ii=jj=1
  for(ii in 1:length(lambda))
  {
    for(jj in 1:length(gamma))
    {
      cYn = Yn - Ymean
      cXn = Xn - matrix(rep(Xmean,dim(Xn)[1]),byrow=TRUE,nrow=dim(Xn)[1])

      slosresult = slos_temp(Y=cYn, X=cXn, Maxiter=Maxiter,lambda=lambda[ii],gamma=gamma[jj],beta.basis=beta.basis,absTol=absTol,Cutoff=Cutoff)
      beta.temp =  slosresult$beta
      beta.slos.all[,ii,jj] =  eval.fd(knots,beta.temp)
      b.slos.all[,ii,jj] = slosresult$b

      BIC_temp = BIC_fun(Y=Yn,b=slosresult$b,U=slosresult$U,V=V, n=length(Yn), alpha=gamma[jj])
      BIC_slos[ii,jj]  = BIC_temp$BIC2
    }
  }

  idx = which(BIC_slos == min(BIC_slos), arr.ind = TRUE)

  Optlambda = lambda[idx[1]]
  Optgamma  = gamma[idx[2]]

  beta_slos  = beta.slos.all[,idx[1],idx[2]]
  b_slos     = b.slos.all[,idx[1],idx[2]]
  T1 = T1.hat(beta_slos)

  result = list(beta=beta_slos,extra=list(X = X, Y = Y, M = M, d = d, domain = domain,b = b_slos, delta=T1,Optgamma=Optgamma,Optlambda=Optlambda))
  class(result) = 'slos'
  result
}



# *** SLoS method (Lin et al 2017)

############################################
# *** SLoS Method
############################################
# y: vector, length n, centered response
# x: matrix, n by p, centered
# Maxiter: a numeric number, the maximum number of iterations for convergence of beta
# lambda: a positive numeric number, tuning parameter for fSCAD penalty
# gamma: a positive numeric number, tuning parameter for the roughness penalty
# beta.basis: basis for beta(t), the coefficient function
# absTol: a numeric number, of max(norm(bHat)) is smaller than absTol, we stop another iteration
# Cutoff: a numeric number, if bHat is smaller than Cutoff, set it to zero to avoide numerical unstable
slos_temp = function(Y, X, Maxiter,lambda,gamma,beta.basis,absTol,Cutoff)
{
  n = length(Y)
  m = dim(X)[2]

  # beta b-spline basis
  rng = getbasisrange(beta.basis)
  breaks = c(rng[1],beta.basis$params,rng[2])
  L = beta.basis$nbasis
  M = length(breaks) - 1
  norder = L-M+1
  d = L-M

  L2NNer = sqrt(M/(rng[2]-rng[1]))

  # calculate design matrix U and roughness penalty matrix V
  U = GetDesignMatrix(X=X,beta.basis=beta.basis)
  V = eval.penalty(beta.basis,int2Lfd(2))
  VV = n*gamma*V

  # calculate W
  W = slos.compute.weights(beta.basis)

  # initial estimate of b
  bHat = solve(t(U)%*%U+VV)%*%t(U)%*%Y

  
  if(lambda > 0)
  {
    changeThres = absTol
    bTilde = slosLQA(U,Y,V,bHat,W,gamma,lambda,Maxiter,M,L,L2NNer,absTol,a=3.7)
    bZero = (abs(bTilde) < Cutoff)
    bTilde[bZero] = 0

    bNonZero = !bZero

    U1 = U[,bNonZero]
    V1 = VV[bNonZero,bNonZero]
    bb = solve(t(U1)%*%U1+V1,t(U1)%*%Y)
    bTilde = matrix(0,dim(U)[2],1)
    bTilde[bNonZero,1] = matrix(bb,length(bb),1)
  }

  bNonZero = as.vector((bTilde != 0))

  projMat = U1 %*% solve(t(U1)%*%U1+V1,t(U1))
  result = list(beta=NULL,projMat=projMat,intercept=0,fitted.values=projMat%*%Y)

  betaobj = list()
  bfd = fd(coef=bTilde,basisobj=beta.basis)
  betaobj = bfd

  result$b = bTilde
  result$beta = betaobj
  result$U = U
  class(result) = 'slos'
  result
}

GetDesignMatrix = function(X,beta.basis)
{
  rng = getbasisrange(beta.basis)
  breaks = c(rng[1],beta.basis$params,rng[2])
  M = length(breaks) - 1

  beta.basismat = eval.basis(breaks,beta.basis)
  hDM = (rng[2]-rng[1])/M
  cefDM = c(1, rep(c(4,2), (M-2)/2), 4, 1) # https://en.wikipedia.org/wiki/Simpson%27s_rule
  U = hDM/3*X%*%diag(cefDM)%*%beta.basismat
  return(U=U)
}


slosLQA = function(U,Y,V,bHat,W,gamma,lambda,Maxiter,M,L,L2NNer,absTol,a)
{
  d=3
  betaNormj = c(0,M)
  bZeroMat = rep(FALSE,L)
  betaNorm = Inf
  n = length(Y)

  it = 1

  while(it <= Maxiter)
  {
    betaNormOld = betaNorm
    betaNorm = sqrt(sum(bHat^2))

    change = (betaNormOld-betaNorm)^2
    if(change < absTol) break

    lqaW = NULL
    lqaWk = matrix(0,L,L)
    for(j in 1:M)
    {
      index = j:(j+d)
      betaNormj[j] = t(bHat[j:(j+d)])%*%W[,,j]%*%bHat[j:(j+d)]
      cjk = Dpfunc(betaNormj[j]*L2NNer,lambda,a)

      if(cjk != 0)
      {
        if(betaNormj[j] < absTol) {
          bZeroMat[index] = TRUE
          }else{lqaWk[index,index] = lqaWk[index,index] + cjk*(L2NNer/betaNormj[j])*W[,,j]}
      }
    }

    lqaW = lqaWk
    lqaW = lqaW / 2

    bZeroVec = bZeroMat
    bNonZeroVec = !bZeroVec

    UtU = t(U[,bNonZeroVec])%*%U[,bNonZeroVec]
    Ut = t(U[,bNonZeroVec])
    Vp = n*gamma*V[bNonZeroVec,bNonZeroVec]

    theta = solve(UtU+Vp+n*lqaW[bNonZeroVec,bNonZeroVec,drop=F],Ut %*% Y)
    bHat = matrix(0,length(bNonZeroVec),1)
    bHat[bNonZeroVec] = theta

    it = it + 1
  }
  bHat
}


slos.compute.weights = function(basis)
{
  L = basis$nbasis
  rng = getbasisrange(basis)
  breaks = c(rng[1],basis$params,rng[2])
  M = length(breaks) - 1
  norder = L-M+1
  W = array(0,dim=c(norder,norder,M))
  for (j in 1:M)
  {
    temp = inprod(basis,basis,rng=c(breaks[j],breaks[j+1]))
    W[,,j] = temp[j:(j+norder-1),j:(j+norder-1)]
  }
  W
}

Dpfunc = function(u,lambda,a)
{
  if(u<=lambda) Dpval = lambda
  else if(u<a*lambda) Dpval = -(u-a*lambda)/(a-1)
  else Dpval = 0
  Dpval
}


BIC_fun = function(Y, b, U, V, n, alpha)
{
  sparse.idx   = which(b == 0)
  if(length(sparse.idx) == 0)
  {
    ula = U
    vla = V
  }
  else{
    ula  = U[, -sparse.idx]
    vla  = V[-sparse.idx, -sparse.idx]
  }

  hat2 = ula%*%solve(t(ula)%*%ula + n*alpha*vla)%*%t(ula)
  df2  = tr(hat2)
  Yhat = U%*%b
  RSS.gbric  = t(Y - Yhat)%*%(Y - Yhat)
  BIC2temp.gbr = n*log(RSS.gbric/n) + log(n)*df2
  return(list(BIC2 = BIC2temp.gbr))
}



