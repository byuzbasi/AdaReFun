
#' Nested Group Bridge Regression Model
#'
#' Calculates a functional regression model using a nested group bridge approach.
#'
#' @import fda
#' @import psych
#' @import glmnet
#' @import stats
#' @import utils
#' @import MASS
#' @import flare
#' @param Y Vector of length n.
#' @param X Matrix of n x p, covariate matrix, should be dense.
#' @param M Integer, t1,..., tM are M equally spaced knots.
#' @param d Integer, the degree of B-Splines.
#' @param domain The range over which the function X(t) is evaluated and the coefficient function \eqn{\beta}(t) is expanded by the B-spline basis functions.
#' @param extra List containing other parameters which have defaults:
#'              \itemize{
#'                \item alphaPs: Smoothing parameter for the Penalized B-splines method, default is 10^(-10:0).
#'                \item kappa: Tuning parameter for roughness penalty, default is 10^(-(9:7)).
#'                \item tau: Tuning parameter for the group bridge penalty, default is exp(seq(-50,-15,len = 20)).
#'                \item gamma: Real number, default is 0.5.
#'                \item niter: Integer, maximum number of iterations, default is 100.
#'              }
#'
#' @return beta: Estimated \eqn{\beta}(t) at discrete points.
#' @return extra: List containing other values which may be of use:
#'        \itemize{
#'          \item b: Estimated b-hat.
#'          \item delta: Estimated cutoff point.
#'          \item Ymean: Estimated y-hat.
#'          \item Xmean: Estimated x-hat.
#'          \item Optkappa: Optimal roughness penalty selected.
#'          \item Opttau: Optimal group bridge penalty selected.
#'          \item M: Integer representing the number of knots used in the model calculation.
#'          \item d: Integer, degree of B-Splines used.
#'          \item domain: The range over which the function X(t) was evaluated and the coefficient function \eqn{\beta}(t) was expanded by the B-spline basis functions.
#'        }
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
#'
ngr = function(Y, X, M,d,domain, extra=list(alphaPS=10^(-10:0), kappa=10^(-(9:7)), tau=exp(seq(-50,-15,len = 20)), gamma=0.5,niter=100))
{
  p = dim(X)[2]
  norder   = d+1
  knots    = seq(domain[1],domain[2], length.out=M+1)
  nknots   = length(knots)
  nbasis   = length(knots) + norder - 2 # i+2
  basis  = create.bspline.basis(knots,nbasis,norder)
  V = eval.penalty(basis,int2Lfd(2))
  tobs = seq(domain[1],domain[2],length.out = p)
  basismat = eval.basis(tobs, basis)

  if (is.null(extra$alphaPS)) {
    alphaPS=exp(seq(-20,20,len = 20))
  } else {
    alphaPS=extra$alphaPS
  }

  if (is.null(extra$kappa)) {
    kappa=exp(seq(-10,10,len = 8))
  } else {
    kappa=extra$kappa
  }

  if (is.null(extra$tau)) {
    tau=exp(seq(-50,-10,len = 8))
  } else {
    tau=extra$tau
  }

  if (is.null(extra$gamma)) {
    gamma=0.5
  } else {
    gamma=extra$gamma
  }

  if (is.null(extra$niter)) {
    niter=100
  } else {
    niter=extra$niter
  }

  n = length(Y)

  # use BIC tune
  nkappa = length(kappa)
  ntau   = length(tau)
  bic.ngr = array(NA,c(nkappa, ntau))
  for(i in 1:nkappa)
  {
    for(j in 1:ntau)
    {
      ngr.fit = ngr.algorithm(Y=Y, X=X, V=V, domain=domain, basismat= basismat, alphaPS=alphaPS, kappa=kappa[i], tau=tau[j], gamma=gamma, M=M, d=d, niter=niter)
      Yhatngr = ngr.fit$fitted.values
      bic.ngr[i,j] = AICBIC.ngr(Y=Y, b=ngr.fit$b, Yhat=Yhatngr, U=ngr.fit$U, V=V, n=n, kappa=kappa[i])$BIC
    }
  }
  idx = which(bic.ngr == min(bic.ngr), arr.ind = TRUE)
  kappa.ngr = kappa[idx[1]]
  tau.ngr = tau[idx[2]]

  fit = ngr.algorithm(Y=Y, X=X, V=V, domain=domain,basismat= basismat, alphaPS=alphaPS, kappa=kappa.ngr, tau=tau.ngr, gamma=gamma, M=M, d=d, niter=niter)


  result = list(beta=fit$beta,extra=list(b=fit$b,delta=fit$delta,Ymean=fit$Ymean,Xmean=fit$Xmean,Optkappa=kappa.ngr,Opttau=tau.ngr,M=M,d=d,domain=domain))
  class(result) = 'ngr'
  result
}









# *** nested group bridge method
# *** penalized B-splines method (Cardot et al 2003)

#########################################################################
# *** Nested Group Bridge Approach + Penalized B-splines Method
#########################################################################
# INPUTs
#   Y: vector of length n
#   X: matrix of n x p, covariate matrix,should be dense
#   U: matrix of n x M+d
#   v: matrix of M+d x M+d, roughness matrix
#   M: integer, t1,...tM are n equally spaced knots
#   d: integer, the degree of Bsplines
#   kappa: real number, tuning parameter for roughness penalty
#   tau: real number, tuning parameter for the group bridge penalty
#   gamma: real number, use 0.5
Ystar = function(Y, alpha, M, d)
{if(alpha == 0){tempystar = Y}else{tempystar = as.matrix(rbind(as.matrix(Y, ncol = 1), matrix(0, nrow = M + d, ncol = 1)))}
  return(tempystar)
}

Ustar = function(U, V, n, alpha, M, d)
{
  if(alpha==0){tempustar = U}else{eig = eigen(V)
  eig$values[eig$values < 0] = 0
  W   = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)

  tempustar = rbind(U, sqrt(n*alpha)*W)}
  return(tempustar)
}

# penalized B-splines
PS_ngr = function(Y, U, V, n, alpha, M, d)
{
  Usmoothb  = Ustar(U = U, V = V, n = n, alpha = alpha, M = M, d = d)
  Ysmoothb  = Ystar(Y = Y, alpha = alpha, M = M, d = d)

  smoothfit = lm(Ysmoothb ~  0 + Usmoothb)
  Ystarhat  = smoothfit$fitted.values
  bsmoothhat = smoothfit$coefficients

  return(list(bsmoothhat = bsmoothhat, Ysmoothhat = Ystarhat[1:n]))
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

# alpha: a vector
PS.tune.BIC = function(Y, U, V, n, alpha, M, d, plot.it=TRUE)
{
  # use BIC tune
  nalpha = length(alpha)
  bic.ps = rep(NA,nalpha)
  for(i in 1:nalpha)
  {
    Yhatps    = PS_ngr(Y=Y, U=U, V=V, n = n, alpha = alpha[i], M, d)$Ysmoothhat
    bic.ps[i] = AICBIC.PS(Y=Y, Yhat=Yhatps, U=U, V=V, n=n, alpha=alpha[i])$BIC
  }

  alpha.ps = alpha[which.min(bic.ps)]


  if(plot.it)
  {
    #plot(alpha,bic.ps,type="b",col=2,pch=16,lwd=1.5,xlab="alpha",ylab="BIC")
    #abline(v=alpha.ps,lty=2,col=2);axis(3,alpha.ps)
  }
  return(list(Optalpha=alpha.ps,bic=bic.ps))
}



# Nested group bridge method
theta_js = function(b_sminus1, cj, n, gamma, tau, M, d, j)
{
  b2bargamma = (sum(abs(b_sminus1)[j:(d+M)]))^gamma
  tempthetajs = cj*((1/gamma-1)/tau)^gamma*b2bargamma
  return(tempthetajs)
}

h_js = function(theta_js, cj, gamma)
{
  temphjs = theta_js^(1-1/gamma)*cj^(1/gamma)
  return(temphjs)
}

g_ks = function(h_js, k, M, d)
{
  if(k <= M) {tempgks = sum(h_js[1:k])}
  else {tempgks = sum(h_js[1:M])}
}


cj = function(M, d, gamma, j, bsmooth)
{
  tempcj  = (d + M + 1 - j)^(1-gamma)
  bsm_Aj_norm_gamma  = sqrt(sum((bsmooth[j:(M+d)])^2))^gamma
  tempcjfinal = tempcj/bsm_Aj_norm_gamma
  return(tempcjfinal)
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

compute.u = function(X, n, M, d, domain)
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



deltaind = function(b,p)
{
  ind = 0
  indi = 0
  sumbtemp = rep(NA,length(b))
  indtemp = rep(NA,length(b))
  for(i in 1:length(b))
  {
    sumbtemp[i] = sum(abs(b[i:length(b)]))
    if(sumbtemp[i]==0){indtemp[i]=1}else{indtemp[i]=0}
  }
  if(sum(indtemp)==0) {indi=p}else{indi=which(indtemp==1)[1]}
  return(indi)
}


ngr.algorithm = function(Y, X, V, domain, basismat, alphaPS=exp(seq(-20,20,len = 20)), kappa=exp(seq(-10,10,len = 8)), tau=exp(seq(-50,-10,len = 8)), gamma=0.5, M=100, d=3, niter=100)
{
  tobs = seq(domain[1], domain[2], length.out = dim(X)[2])
  n=length(Y)
  p = dim(X)[2]
  yc = Y.center(Y)$Ycenter
  xc = X.center(X)$Xcenter
  uc = compute.u(X=xc, n=n, M=M, d=d, domain=domain)$U

  pstune = PS.tune.BIC(Y=yc, U=uc, V=V, n=n, alpha=alphaPS, M=M, d=d, plot.it=FALSE)
  psfit = PS_ngr(Y=yc, U=uc, V=V, n=n, alpha=pstune$Optalpha, M=M, d=d)
  b0 = psfit$bsmoothhat
  #betaPS =  basismat%*%b0

  biter = matrix(NA, nrow = length(b0), ncol = niter)
  for(iterbg in 1:niter)
  {
    if(iterbg == 1) {biter[,iterbg] = b0}
    ##################################
    #    Step 1: compute theta_js    #
    ##################################
    else {
      theta = rep(0, M)
      hn     = rep(0, M)
      for(j in 1:M)
      {
        theta[j] = theta_js(b_sminus1 = biter[,(iterbg-1)], cj = cj(M = M, d = d, gamma = gamma, j = j, bsmooth = b0), n = n, gamma = gamma, tau = tau, M = M, d = d, j = j)
        hn[j]     = h_js(theta[j], cj(M = M, d = d, gamma = gamma, j = j, bsmooth = b0), gamma)
      }


      ##################################
      #    Step 2: compute g_ks        #
      ##################################
      g = rep(0, (M + d))
      for(k in 1:(M + d))
      {
        g[k] = g_ks(h_js = hn, k = k, M = M, d = d)
      }

      ##################################
      #    Step 3: compute bs          #
      ##################################
      Ustarbhatgbr = Ustar(U = uc, V = V, n = n, alpha = kappa, M = M, d = d)
      Ystarbhatgbr = Ystar(Y = yc, alpha = kappa, M = M, d = d)
      Ustarstargbr = Ustarbhatgbr%*%diag(1/(n*g))
      lassomodel   = glmnet(x = Ustarstargbr, y = Ystarbhatgbr, standardize = FALSE, alpha = 1, lambda = 0.5/n, family = "gaussian", intercept = FALSE)
      biter[,iterbg] =  coef(lassomodel)[2:length(coef(lassomodel)),]/(n*g)

      ################################################
      # decide whether to stop the iteration         #
      ################################################

      difratio = rep(0, length(biter[,iterbg]))
      if(iterbg >= 3){
        idx0 = which((biter[,iterbg]-biter[,(iterbg-1)]) == 0)
        if(length(idx0) == 0){difratio = (biter[,iterbg] - biter[, (iterbg-1)])/biter[,(iterbg-1)]}
        else {difratio[-idx0] = (biter[-idx0,iterbg] - biter[-idx0, (iterbg-1)])/biter[-idx0,(iterbg-1)]}
        if(max(difratio) < 10^(-4)) break}
    }
  }


  bhatgbr = biter[,iterbg]
  Yhat = uc%*%bhatgbr + mean(Y)
  deltaindicator = deltaind(bhatgbr,p)
  deltangr = tobs[deltaindicator]

  return(list(beta=basismat%*%bhatgbr,b=bhatgbr,delta=deltangr,fitted.values=Yhat,U=uc,Ymean=mean(Y),Xmean=apply(X,2,mean)))
}


AICBIC.ngr = function(Y, b, Yhat, U, V, n, kappa)
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

  hat = ula%*%solve(t(ula)%*%ula + n*kappa*vla)%*%t(ula)
  df  = tr(hat)

  SSE = t(Y - Yhat)%*%(Y - Yhat)
  AIC.ngr = n*log(SSE/n) + 2*df
  BIC.ngr = n*log(SSE/n) + log(n)*df
  return(list(AIC = AIC.ngr, BIC = BIC.ngr))
}












beta_fun = function(t, ii)
{
  if(ii == 1) {bf = ifelse(t<=0.5 && t>=0,1,0)}
  else if(ii == 2){bf = sin(2*pi*t)*ifelse(t<=0.5 && t>=0,1,0)}
  else if(ii == 3){bf = (cos(2*pi*t)+1)*ifelse(t<=0.5 && t>=0,1,0)}
  else if(ii == 4){bf = (-100*(t-0.5)^3 - 200*(t - 0.5)^4)*ifelse(t<=0.5 && t>=0,1,0)}
  else{print("model does not exit")}
  return(bf)
}



#' Generating random curves from B-Splines
#'
#' Generating random curves from B-Splines
#'n,nknots,norder,p,domain=c(0,1),snr,betaind
#' @param n Number of curves
#' @param nknots Number of knots
#' @param norder Degree
#' @param p Number of time points
#' @param domain Domain of time
#' @param snr Signal to noise ratio
#' @param betaind Numeric index for function
#'
#' @return X: The generated X matrix of curve sampled at each timepoint
#' @return Y: The generated dependent variable
#' @export
ngr.data.generator.bsplines = function(n,nknots,norder,p,domain=c(0,1),snr,betaind)
{
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2
  basis    = create.bspline.basis(knots,nbasis,norder)

  tobs = seq(domain[1],domain[2],length.out = p)
  basismat = eval.basis(tobs, basis)

  x=array(NA,c(n,p))
  for(i in 1:n)
  {
    x[i,] = rnorm(nbasis, 0, 1)%*%t(basismat)
  }

  betaeval = apply(as.matrix(tobs), 1, beta_fun, ii = betaind)

  # y0 the signals
  h   = (domain[2]-domain[1])/(p-1)
  cef = c(1, rep(c(4,2), (p-3)/2), 4, 1)
  y0  = rep(NA,n)
  y0  = h/3*x%*%diag(cef)%*%betaeval
  eps0= sd(y0)
  y   = y0 + rnorm(n,mean = 0, sd = eps0)

  return(list(X=x,Y=y))
}




ngr.data.generator.fourier = function(n,nfourier.basis,p,var.alpha,domain=c(0,1),snr,betaind)
{
  tobs = seq(domain[1],domain[2],length.out = p)
  fbasis = create.fourier.basis(domain,nbasis=nfourier.basis)
  psi = eval.basis(tobs,fbasis)

  x =  psi%*%diag(1/(1:nfourier.basis)^var.alpha)%*%matrix(rnorm(nfourier.basis*n),nfourier.basis,n) # relatively rough

  betaeval = apply(as.matrix(tobs), 1, beta_fun, ii = betaind)

  # y0 the signals
  h   = (domain[2]-domain[1])/(p-1)
  cef = c(1, rep(c(4,2), (p-3)/2), 4, 1)
  y0  = rep(NA,n)
  y0  = h/3*t(x)%*%diag(cef)%*%betaeval
  eps0= sd(y0)
  y   = y0 + rnorm(n,mean = 0, sd = eps0)

  return(list(X=t(x),Y=y))
}


# define function for the estimator given alpha, tau, gamma and niter
# INPUT: 1, Y: response, vector of length n
#        2, alpha : tunning parameter for roughness penalty
#        2, lambda : tunning parameter ralated to group bridge penalty term with gamma = 1
#        3, M: t0,t1,...tM are konts
#        4, d: degree of freedom of B-spline basis
#        5. u: n x M+d matrix
#        6. v: M+d x M+d matrix
#        7. b0: used to calculate cj


# OUTPUT: 1, bhat: estimator
#         2, RSS : RSS of the validation sets
#         3. betahat: coefficient function
# grouop bridge AIC and BIC
AICBIC.gbr.fun2 = function(Y, b, U, V, n, alpha)
{
  sparse.idx   = which(b == 0)
  if(length(sparse.idx) == 0)
  {
    ula = U
    vla = V
    #wla = W
  }
  else{
    ula  = U[, -sparse.idx]
    vla  = V[-sparse.idx, -sparse.idx]
    #wla  = W[-sparse.idx, -sparse.idx]
  }
  #hat1 = ula%*%solve(t(ula)%*%ula + n*alpha*vla + 0.5*wla)%*%t(ula)
  hat2 = ula%*%solve(t(ula)%*%ula + n*alpha*vla)%*%t(ula)
  #df1  = tr(hat1)
  df2  = tr(hat2)
  Yhat = U%*%b
  RSS.gbric  = t(Y - Yhat)%*%(Y - Yhat)
  #AIC1temp.gbr = n*log(RSS.gbric/n) + 2*df1
  AIC2temp.gbr = n*log(RSS.gbric/n) + 2*df2
  #BIC1temp.gbr = n*log(RSS.gbric/n) + log(n)*df1
  BIC2temp.gbr = n*log(RSS.gbric/n) + log(n)*df2
  return(list(AIC2 = AIC2temp.gbr, BIC2 = BIC2temp.gbr))
}

weightedlasso = function(Y, X, M,d,domain,V, n, alpha,lambda,alphaPS,basismat)
{
  p = dim(X)[2]
  yc = Y.center(Y)$Ycenter
  xc = X.center(X)$Xcenter
  uc = compute.u(X=xc, n=n, M=M, d=d, domain=domain)$U

  pstune = PS.tune.BIC(Y=yc, U=uc, V=V, n=n, alpha=alphaPS, M=M, d=d, plot.it=FALSE)
  psfit = PS_ngr(Y=yc, U=uc, V=V, n=n, alpha=pstune$Optalpha, M=M, d=d)
  b0 = psfit$bsmoothhat

  cj=rep(NA,M+d)
  for(iterc in 1:(M+d))
  {
    cj[iterc] = 1/sqrt(sum((b0[iterc:(M+d)])^2))
  }

  ##################################
  #    Step 1: compute g_ks        #
  ##################################
  g = rep(0, (M + d))
  for(k in 1:(M + d))
  {
    g[k] = sum(cj[1:min(c(k,M))])
  }

  ##################################
  #    Step 3: compute bs          #
  ##################################
  Ustarbhatgbr = Ustar(U = uc, V = V, n = n, alpha = alpha, M = M, d = d)
  Ystarbhatgbr = Ystar(Y = Y, alpha = alpha, M = M, d = d)
  Ustarstargbr = Ustarbhatgbr%*%diag(1/(lambda*n*g))
  lassomodel   = glmnet(x = Ustarstargbr, y = Ystarbhatgbr, standardize = FALSE, alpha = 1, lambda = 0.5/n, family = "gaussian", intercept = FALSE)
  bhatgbr =  coef(lassomodel)[2:length(coef(lassomodel)),]/(lambda*n*g)
  tobs=seq(domain[1],domain[2],length.out = p)
  Yhat = uc%*%bhatgbr + mean(Y)
  deltaindicator = deltaind(bhatgbr,p)
  deltangr = tobs[deltaindicator]

  result = list(b.wl=bhatgbr,beta.wl=basismat%*%bhatgbr,delta.wl=deltangr,fitted.values=Yhat, Ymean=mean(Y),Xmean=apply(X,2,mean),U=uc, bPS=b0)
  class(result) = 'wl'
  result
}


weightedlasso.tune = function(Y, X, V, n, alphaPS, alpha,lambda, M, d,domain,basismat)
{
  nalpha = length(alpha)
  nlambda   = length(lambda)
  bic.wl = array(NA,c(nalpha, nlambda))
  for(i in 1:nalpha)
  {
    for(j in 1:nlambda)
    {
      wl.fit = weightedlasso(Y=Y,X=X,M=M,d=d,domain=domain,V=V,n=n,alpha=alpha[i],lambda=lambda[j],alphaPS,basismat)
      Yhatwl = wl.fit$fitted.values
      bic.wl[i,j] = AICBIC.gbr.fun2(Y=Y,b=wl.fit$b.wl,U=wl.fit$U,V=V,n=n,alpha=alpha[i])$BIC
    }
  }
  idx = which(bic.wl == min(bic.wl), arr.ind = TRUE)
  alpha.ngr = alpha[idx[1]]
  lambda.ngr = lambda[idx[2]]
  return(list(Optkappa=alpha.ngr,Opttau=lambda.ngr,bic=bic.wl))
}



T1.hat = function(beta)
{
  T1n = length(beta)
  t1hatfuntemp=1
  if(sum(abs(beta)[1:T1n])==0){t1hatfuntemp=0}else{
    for(t1i in 1:(T1n-1))
    {
      if(beta[t1i]!=0 && sum(abs(beta)[(t1i+1):T1n])==0)
      {t1hatfuntemp = t1i/(T1n-1)}
    }}
  return(t1hatfuntemp)
}


