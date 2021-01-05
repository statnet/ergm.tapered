# Kurtosis related functions
#' @export
eq.fun.kurtosis <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){

  Var.xsim <- grep("Var(",colnames(xsim), fixed=TRUE)
# Standardize the Var and M4 against original
  Kurt.xsim.stats <- substr(colnames(xsim)[Var.xsim],5,nchar(colnames(xsim)[Var.xsim])-1)
  Kurt.xsim.stats <- match(Kurt.xsim.stats,colnames(xsim))
  xsim[,Var.xsim] <- xsim[,Kurt.xsim.stats]^2
# xsim[,Var.xsim] <- sweep(xsim[,Kurt.xsim.stats]^2,2,colMeans(xsim[,Kurt.xsim.stats]^2),"-")

  Ki <- control.llik$MCMLE.kurtosis.location

  eta <- ergm.eta(theta, etamap)
  eta[is.na(eta)] <- 0

  etaparam <- (eta-eta0)
  basepred <- xsim %*% etaparam
  maxbase <- max(basepred)
  lwi <- basepred - (maxbase + log(sum(exp(basepred-maxbase))) )
  Ek <- colSums(sweep(xsim[,Kurt.xsim.stats],1,exp(lwi),"*"))
## vardev <- sweep(xsim[,Kurt.xsim.stats],2,Ek,"-")^2
## varj <- colSums(sweep(vardev,1,exp(lwi),"*"))
#  varj <- colSums(sweep(xsim[,Var.xsim],1,exp(lwi),"*"))
## a <- eta[M4.xsim]*(2*Ki*varj) - (-eta[Var.xsim])
## a <- eta[M4.xsim]*(4*Ki*(varj^(1/3))) - (-3*eta[Var.xsim])
## a[is.nan(a)] <- eta[Var.xsim][is.nan(a)]
  logm2 <- log(colSums(sweep(xsim[,Var.xsim]  ,1,exp(lwi),"*")))
  logm4 <- log(colSums(sweep(xsim[,Var.xsim]^2,1,exp(lwi),"*")))
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
  logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
  kurt <- exp(logm4-2*logm2)
# kurt <- logm4-2*logm2
  kurt[is.nan(kurt) | is.na(kurt)] <- 3
  names(kurt) <- colnames(xsim)[Var.xsim]
#print(mean(kurt))
  penalty <- -0.5*((kurt-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
  penalty <- sum(penalty)
# cat(sprintf("mean(kurt)=%f, penalty = %9.5g\n",mean(kurt),penalty))
# print(round(kurt,2))
# c(Ek,kurt-log(Ki))
  c(Ek,kurt-Ki)
}
#end

#' @export
llik.fun.kurtosis <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){
# Convert theta to eta
  Var.xsim <- grep("Var(",colnames(xsim), fixed=TRUE)
# Standardize the Var and M4 against original
  Kurt.xsim.stats <- substr(colnames(xsim)[Var.xsim],5,nchar(colnames(xsim)[Var.xsim])-1)
  Kurt.xsim.stats <- match(Kurt.xsim.stats,colnames(xsim))
  xsim[,Var.xsim] <- xsim[,Kurt.xsim.stats]^2
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  basepred <- xsim %*% etaparam
  mb <- lweighted.mean(basepred,rowweights(xsim))
  vb <- lweighted.var(basepred,rowweights(xsim))
# This is log-lik (and not its negative)
  llr <- -mb - control.llik$MCMLE.varweight*vb
#if(llr < -800){browser()}
#
  penalty <- eta[Var.xsim]
  penalty[penalty > 0] <- -10/control.llik$MCMLE.kurtosis.scale - penalty[penalty > 0]
# penalty <- mean(penalty*penalty)
# cat(sprintf("llr = %f, penalty = %9.5g\n",llr, penalty))
# llr <- llr + control.llik$MCMLE.kurtosis.scale*penalty

  # Simplistic error control;  -800 is effectively like -Inf:
  if(is.infinite(llr) | is.na(llr)){llr <- -200}

  maxbase <- max(basepred)
  lwi <- basepred - (maxbase + log(sum(exp(basepred-maxbase))) )
  logm2 <- log(colSums(sweep(xsim[,Var.xsim]  ,1,exp(lwi),"*")))
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
  penalty <- mean(penalty*exp(logm2),na.rm=TRUE)/mean(exp(logm2))
#
# llr <- llr + control.llik$MCMLE.kurtosis.scale*penalty
#
# logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
# kurt <- exp(logm4-2*logm2)
# kurt[is.nan(kurt) | is.na(kurt)] <- 3
# names(kurt) <- colnames(xsim)[Var.xsim]
# penalty <- -0.5*((kurt[-Var.xsim]-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
# penalty <- sum(penalty)
  # trustregion is the maximum value of llr that we actually trust.
  # So if llr>trustregion, return a value less than trustregion instead.
  #
  # Return the negative log-lik (for nloptr is a minimizer)
# if (is.numeric(control.llik$MCMLE.trustregion) && llr< -1000) {
  if (llr < -800) {
    return(800)
  } else {
    return(-llr)
  }
}
#' @export
llik.grad.kurtosis <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){

  Var.xsim <- grep("Var(",colnames(xsim), fixed=TRUE)
# Standardize the Var and M4 against original
  Kurt.xsim.stats <- substr(colnames(xsim)[Var.xsim],5,nchar(colnames(xsim)[Var.xsim])-1)
  Kurt.xsim.stats <- match(Kurt.xsim.stats,colnames(xsim))
  xsim[,Var.xsim] <- xsim[,Kurt.xsim.stats]^2
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights
  basepred <- xsim %*% etaparam + lrowweights(xsim)
  
  # Calculate the estimating function values sans offset
  llg <- - lweighted.mean(xsim, basepred)
  llg <- t(ergm.etagradmult(theta, llg, etamap))
  
  llg[is.na(llg)] <- 0

  # Return negative grad
  -llg
}

#' @export
eq.jacobian.kurtosis <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){

  Var.xsim <- grep("Var(",colnames(xsim), fixed=TRUE)
  npar <- length(Var.xsim)
# Standardize the Var and M4 against original
  Kurt.xsim.stats <- substr(colnames(xsim)[Var.xsim],5,nchar(colnames(xsim)[Var.xsim])-1)
  Kurt.xsim.stats <- match(Kurt.xsim.stats,colnames(xsim))
  xsim[,Var.xsim] <- xsim[,Kurt.xsim.stats]^2

  Ki <- control.llik$MCMLE.kurtosis.location

  eta <- ergm.eta(theta, etamap)
  eta[is.na(eta)] <- 0

  etaparam <- (eta-eta0)
  basepred <- xsim %*% etaparam
  maxbase <- max(basepred)
  lwi <- basepred - (maxbase + log(sum(exp(basepred-maxbase))) )
#
  Ek <- colSums(sweep(xsim,1,exp(lwi),"*"))
  t2 <- sweep(sweep(xsim,2,Ek,"-"),1,exp(lwi),"*")
  jac_Var <- t(xsim[,Var.xsim]) %*% t2
  jac_M4  <- t(xsim[,Var.xsim]^2) %*% t2
  jac0   <-  t(xsim[,Kurt.xsim.stats]) %*% t2
  jac <- rbind(jac0, jac_M4-2*jac_Var)
  jac[is.na(jac) | is.nan(jac)] <- 0
  jac
}

#' @export
llik.hessian.kurtosis <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){

  Var.xsim <- grep("Var(",colnames(xsim), fixed=TRUE)
  M4.xsim <-  grep("M4(", colnames(xsim), fixed=TRUE)
  Kurt.xsim <- c(Var.xsim,M4.xsim)
  eta0 <- eta0[seq_along(theta)]

  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0
  etaparam[Kurt.xsim] <- 0

  # Calculate log-importance-weights
  basepred <- xsim %*% etaparam + lrowweights(xsim)

  # Calculate the estimating function values sans offset
  esim <- t(ergm.etagradmult(theta, t(xsim), etamap))
  esim[,Kurt.xsim] <- 0

  # Weighted variance-covariance matrix of estimating functions ~ -Hessian
  H <- -lweighted.var(esim, basepred)

  dimnames(H) <- list(names(theta), names(theta))
  H
}

#' @export
llik.hessian.kurtosis.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){
H <-  numDeriv::hessian(llik.fun.kurtosis,theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                   )
# cnames <- names(theta)
# M4.names <- cnames[grep("offset(M4(",cnames,fixed=TRUE)]
# var.names <- sapply(strsplit(substr(M4.names,start=11,stop=1000),")",fixed=TRUE),
#                     function(x){x[1]})
# cnames <- match(var.names,cnames)
# H[-cnames,] <- 0
# H[,-cnames] <- 0
# diag(H)[-cnames] <- 1
#
# Note we return the Hessian of the log-lik (and not the negatie log-like, which is that optimized)
-H
}

#' @export
llik.hessian.kurtosis.direct <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){
  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights
  basepred <- xsim %*% etaparam + lrowweights(xsim)

  # Calculate the estimating function values sans offset
  esim <- t(ergm.etagradmult(theta, t(xsim), etamap))

  # Weighted variance-covariance matrix of estimating functions ~ -Hessian
  H <- -lweighted.var(esim, basepred)

  dimnames(H) <- list(names(theta), names(theta))

  kfn <- function(theta, xsim, etamap, eta0){
    # Convert theta to eta
    if(is.null(theta)){
      logm2 <- log(apply(sweep(xsim^2,1,rowweights(xsim),"*"),2,sum))
      logm4 <- log(apply(sweep(xsim^(8/3),1,rowweights(xsim),"*"),2,sum))
      logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
      logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
    }else{
      eta <- ergm.eta(theta, etamap)
      # Calculate approximation to the penalized l(eta) - l(eta0) using a lognormal approximation
      etaparam <- eta-eta0
      basepred <- xsim %*% etaparam
      maxbase <- max(basepred)
      logm2 <- maxbase + log(apply(sweep(xsim^2,1,rowweights(xsim)*exp(basepred-maxbase),"*"),2,sum))
      logm4 <- maxbase + log(apply(sweep(xsim^(8/3),1,rowweights(xsim)*exp(basepred-maxbase),"*"),2,sum))
      logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
      logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
    }
    kurt <- exp(logm4-2*logm2)
    kurt[is.nan(kurt) | is.na(kurt)] <- 3
    exterms <- grep("Var(",colnames(xsim),fixed=TRUE)
    names(kurt) <- colnames(xsim)[-exterms]
    penalty <- -0.5*((kurt[-exterms]-3)/control.llik$MCMLE.kurtosis.scale)^2
    penalty <- sum(penalty)
#   if(verbose) {
#     message("kurtosis=")
#     print(kurt[-exterms])
#   }
    penalty
  }

  if(F){
   hessian.penalty <- numDeriv::hessian(kfn, theta, 
    method.args=list(eps=1e-2, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE),
    xsim=xsim, etamap=etamap, eta0=eta0)
   #print(dim(H))
   #print(dim(hessian.penalty))
   print((H))
   print((hessian.penalty))
    H+hessian.penalty
  }else{
    H
  }
}

#llik.grad.kurtosis.orig <- llik.grad.IS
#llik.hessian.kurtosis.orig <- llik.hessian.IS
# Numerical versions
#' @export
llik.grad.kurtosis.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){
  numDeriv::grad(llik.fun.kurtosis,theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                 )
}
#' @export
eq.jacobian.kurtosis.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){
# numDeriv::jacobian(eq.fun.kurtosis,theta,
  nloptr::nl.jacobian(x0=theta,fn=eq.fun.kurtosis,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                 )
}
