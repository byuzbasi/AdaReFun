
#' Penalized B-splines Regression Model
#'
#' Calculates a functional regression model using the penalized B-splines method.
#'
#' @import fda
#' @import psych
#' @import glmnet
#' @import stats
#' @import utils
#' @import MASS
#' @import flare
#' @param Y Vector of length n.
#' @param X  Matrix of n x p, covariate matrix, should be dense.
#' @param alpha Vector.
#' @param M Integer, t1,..., tM are M equally spaced knots.
#' @param d Integer, the degree of B-Splines.
#' @param domain The range over which the function X(t) is evaluated and the coefficient function \eqn{\beta}(t) is expanded by the B-spline basis functions.
#'
#' @return beta: Estimated \eqn{\beta}(t) at discrete points.
#' @return extra: List containing other values which may be of use:
#'                \itemize{
#'                  \item b: Estimated B-spline coefficients.
#'                  \item Ymean: Mean of the Y values.
#'                  \item Xmean: Mean of all X values.
#'                  \item Optalpha: Optimal alpha value chosen.
#'                  \item M: Integer representing the number of knots used in the model calculation.
#'                  \item d: Integer, degree of B-Splines used.
#'                  \item domain: The range over which the function X(t) was evaluated and the coefficient function \eqn{\beta}(t) was expanded by the B-spline basis functions.
#'                }
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
#'
#'
PenS = function(Y, X, alpha, M, d, domain)
{
  # use BIC tune
  X = t(X)
  nalpha = length(alpha)
  bic.ps = rep(NA,nalpha)
  n = length(Y)

  for(i in 1:nalpha)
  {
    ps.fit    = PS.algorithm(Y=Y, X=X, alpha = alpha[i], M, d, domain)
    Yhatps    = ps.fit$extra$Ysmoothhat
    bic.ps[i] = AICBIC.PS(Y=Y, Yhat=Yhatps, U=ps.fit$extra$U, V=ps.fit$extra$V, n=n, alpha=alpha[i])$BIC
  }

  alpha.ps = alpha[which.min(bic.ps)]
  ps.fit   = PS.algorithm(Y=Y, X=X, alpha = alpha.ps, M, d, domain)

  result = list(beta=ps.fit$beta,extra=list(b=ps.fit$extra$bsmoothhat,Ymean=mean(Y),Xmean=apply(X,1,mean),Optalpha=alpha.ps,M=M,d=d,domain=domain))
  class(result) = 'ps'
  result
}




PS.algorithm = function(Y, X, alpha, M, d, domain)
{
  X = t(X)
  p = dim(X)[2]
  norder   = d+1
  knots    = seq(domain[1],domain[2], length.out=M+1)
  nknots   = length(knots)
  nbasis   = length(knots) + norder - 2 # i+2
  basis  = create.bspline.basis(knots,nbasis,norder)
  V = eval.penalty(basis,int2Lfd(2))
  tobs = seq(domain[1],domain[2],length.out = p)
  basismat = eval.basis(tobs, basis)

  n=length(Y)

  yc = Y.center(Y)$Ycenter
  xc = X.center(X)$Xcenter
  U = compute.u_PS(X=xc, M=M, d=d, domain=domain)$U

  Usmoothb  = Ustar_PS(U = U, V = V, n = n, alpha = alpha)
  Ysmoothb  = Ystar_PS(Y = Y, alpha = alpha, M = M, d = d)

  smoothfit = lm(Ysmoothb ~  0 + Usmoothb)
  Ystarhat  = smoothfit$fitted.values
  bsmoothhat = smoothfit$coefficients

  beta = basismat%*%bsmoothhat

  result = list(beta=beta, extra=list(Y=Y, X=X, M=M, d=d, V = V, U = U, domain=domain, bsmoothhat = bsmoothhat, Ysmoothhat = Ystarhat[1:n]))
  class(result) = "ps"
  result
}




Ystar_PS = function(Y, alpha, M, d)
{if(alpha == 0){tempystar = Y}else{tempystar = as.matrix(rbind(as.matrix(Y, ncol = 1), matrix(0, nrow = M + d, ncol = 1)))}
  return(tempystar)
}

Ustar_PS = function(U, V, n, alpha)
{
  if(alpha==0){tempustar = U}else{eig = eigen(V)
  eig$values[eig$values < 0] = 0
  W   = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)

  tempustar = rbind(U, sqrt(n*alpha)*W)}
  return(tempustar)
}

Y.center = function(Y){return(list(Ycenter = Y-mean(Y), Ymean = mean(Y)))}

X.center = function(X)
{
  n=dim(X)[1]
  Xmean = apply(X,2,mean)

  matx = matrix(rep(Xmean,n),nrow=n,byrow=T)

  Xcenter  = X - matx
  return(list(Xcenter = Xcenter, Xmean = Xmean))
}

compute.u_PS = function(X, M, d, domain)
{
  norder   = d+1
  nknots   = M+1
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2
  basis    = create.bspline.basis(knots,nbasis,norder)
  p   = dim(X)[2]
  tobs = seq(domain[1], domain[2], length.out = p)
  basismat = eval.basis(tobs, basis) # p x M+d


  h   = (domain[2]-domain[1])/(p-1)
  cef = c(1, rep(c(4,2), (p-3)/2), 4, 1)
  u   = h/3*X%*%diag(cef)%*%basismat
  return(list(U=u))
}


AICBIC.PS = function(Y, Yhat, U, V, n, alpha)
{
  hat.s  = U%*%solve(t(U)%*%U + n*alpha*V)%*%t(U)
  df.s   = tr(hat.s)
  Yhat.s = Yhat
  RSS.s  = t(Y - Yhat.s)%*%(Y - Yhat.s)
  AICtemp.s = n*log(RSS.s/n) + 2*df.s
  BICtemp.s = n*log(RSS.s/n) + log(n)*df.s
  return(list(AIC = AICtemp.s, BIC = BICtemp.s))
}
