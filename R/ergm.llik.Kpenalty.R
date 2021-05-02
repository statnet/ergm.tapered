# Version working on the tapering model with just a penalty
# for Kurtosis and large tau (or small r > 0)
#' @export
eq.fun.Kpenalty <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){

  Kpenalty <- grep("Taper_Penalty",colnames(xsim), fixed=TRUE)

  eta <- ergm.eta(theta, etamap)
  eta[is.na(eta)] <- 0

  etaparam <- (eta-eta0)
  basepred <- xsim %*% etaparam
  maxbase <- max(basepred)
  lwi <- basepred - (maxbase + log(sum(exp(basepred-maxbase))) )
  Ek <- colSums(sweep(xsim[,-Kpenalty],1,exp(lwi),"*"))
  logm2 <- (colSums(sweep(xsim[,-Kpenalty]^2,1,exp(lwi),"*")))
  logm4 <- (colSums(sweep(xsim[,-Kpenalty]^4,1,exp(lwi),"*")))
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 1
  logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 1
# kurt <- exp(logm4-2*logm2)
# kurt <- logm4-2*logm2
  kurt <- logm4/(logm2*logm2)
  kurt[is.nan(kurt) | is.na(kurt)] <- 3
# names(kurt) <- colnames(xsim)[-Kpenalty]
#print(mean(kurt))
  penalty <- -0.5*((kurt-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
  penalty <- mean(penalty)
# cat(sprintf("mean(kurt)=%f, penalty = %9.5g\n",mean(kurt),penalty))
# a=(c(Ek,kurt-control.llik$MCMLE.kurtosis.location))
# if(length(a)!=20 | any(is.na(a)) | any(is.nan(a))){browser()}
# print(as.vector(c(Ek,kurt-control.llik$MCMLE.kurtosis.location)))
# as.vector(c(Ek,kurt-control.llik$MCMLE.kurtosis.location))
#if(any(abs(Ek)>100)) {print(Ek);browser();print(str(Ek))}
#if(any(abs(Ek)>100)) {browser()}
  c(as.vector(Ek), kurt-control.llik$MCMLE.kurtosis.location)
  c(Ek, kurt-control.llik$MCMLE.kurtosis.location)
  mean(kurt)-control.llik$MCMLE.kurtosis.location
}
#end

#' @export
llik.fun.Kpenalty <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){

  Kpenalty <- grep("Taper_Penalty",colnames(xsim), fixed=TRUE)

# Convert theta to eta
  eta <- ergm.eta(theta, etamap)
# if(eta[Kpenalty] < -10){
#   eta[Kpenalty] <- 0
# }else{
#   eta[Kpenalty] <- - 3*exp(eta[Kpenalty]) / (1 + exp(eta[Kpenalty]))
# }
# if(eta0[Kpenalty] < -10){
#   eta0[Kpenalty] <- 0
# }else{
#   eta0[Kpenalty] <- - exp(eta[Kpenalty])
#   eta0[Kpenalty] <- - 3*exp(eta0[Kpenalty]) / (1 + exp(eta0[Kpenalty]))
# }
# eta0[Kpenalty] <- - exp(eta0[Kpenalty])
  etaparam <- eta-eta0

  basepred <- xsim %*% etaparam
  mb <- lweighted.mean(basepred,rowweights(xsim))
  vb <- lweighted.var(basepred,rowweights(xsim))
# This is log-lik (and not its negative)
  llr <- -mb - control.llik$MCMLE.varweight*vb
#
  # Simplistic error control;  -800 is effectively like -Inf:
  if(is.infinite(llr) | is.na(llr)){llr <- -200}
  if(llr < -200){llr <- -200}

  logm2 <- log(colMeans(xsim[,-Kpenalty]^2))
  logm4 <- log(colMeans(xsim[,-Kpenalty]^4))
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
  logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
  kurt <- exp(logm4-2*logm2)
  kurt[is.nan(kurt) | is.na(kurt)] <- 3
# if(kurt < 2) kurt <- 2
# penalty <- -0.5*((kurt-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
# if(kurt > 0.5*control.llik$MCMLE.kurtosis.location){
    penalty <- -0.5*((kurt-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
# }else{
#   penalty <- 0
# }
  Kurt_penalty0 <- mean(penalty)

  maxbase <- max(basepred)
  lwi <- basepred - (maxbase + log(sum(exp(basepred-maxbase))) )
  logm2 <- log(colSums(sweep(xsim[,-Kpenalty]^2,1,exp(lwi),"*")))
  logm4 <- log(colSums(sweep(xsim[,-Kpenalty]^4,1,exp(lwi),"*")))
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
  logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
  kurt <- exp(logm4-2*logm2)
  kurt[is.nan(kurt) | is.na(kurt)] <- 3
# if(kurt < 2) kurt <- 2
# if(length(control.llik$MCMLE.kurtosis.location)==0){browser()}
# if(kurt > 0.5*control.llik$MCMLE.kurtosis.location){
    penalty <- -0.5*((kurt-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
# }else{
#   penalty <- 0
# }
  Kurt_penalty <- mean(penalty)
#
  Kurt_penalty_diff <- Kurt_penalty - Kurt_penalty0
#
  penalty <- -eta[Kpenalty]^2 + eta0[Kpenalty]^2
#if(eta[Kpenalty]!=eta0[Kpenalty]) browser()
# cat(sprintf("llr=%f, penalty = %9.5g\n",llr,penalty))
  Ek  <- eta[Kpenalty]*sum(xsim[,Kpenalty]*exp(lwi))
# Ek0 <- eta0[Kpenalty]*mean(xsim[,Kpenalty])
  Ek0 <- eta0[Kpenalty]*sum(xsim[,Kpenalty]*exp(lwi))
  penalty <- Ek - Ek0
#
# if(eta0[Kpenalty] > 0){
#  eta0_Kpenalty <- 0
# }else{
#  eta0_Kpenalty <- eta0[Kpenalty]
# }
# if(eta[Kpenalty] > 0){
#  eta_Kpenalty <- 0
# }else{
#  eta_Kpenalty <- eta[Kpenalty]
# }
# Tpenalty_diff <- eta_Kpenalty - eta0_Kpenalty

  Tpenalty_diff <- eta[Kpenalty] - eta0[Kpenalty]

#
# cat(sprintf("llr=%f, Ek = %9.5g, sum = %9.5g, penalty = %9.5g\n",llr,Ek,sum(xsim[,Kpenalty]*exp(lwi)),penalty))
 if(control.llik$MCMLE.dampening.min.ess < 0){
   cat(sprintf("Kurtosis of statistics:\n"))
   print(kurt)
#  cat(sprintf("kurt %f %f \n",kurt[1], kurt[2]))
   cat(sprintf("eta = %9.5g, eta0 = %9.5g, Kurt_penalty = %9.5g, Kurt_penalty0 = %9.5g\n",
               eta[Kpenalty],eta0[Kpenalty],Kurt_penalty,Kurt_penalty0))
   cat(sprintf("MCMLE.kurtosis.location %f MCMLE.kurtosis.scale %f llr=%f, Kurt_penalty_diff = %9.5g, Tpenalty_diff = %9.5g\n",control.llik$MCMLE.kurtosis.location,control.llik$MCMLE.kurtosis.scale,llr,
               Kurt_penalty_diff,Tpenalty_diff))
 }
#
  llr <- llr + control.llik$MCMLE.kurtosis.penalty*Tpenalty_diff + Kurt_penalty_diff
# llr <- llr + Ek
#
# logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
# kurt <- exp(logm4-2*logm2)
# kurt[is.nan(kurt) | is.na(kurt)] <- 3
# names(kurt) <- colnames(xsim)[Var.xsim]
# penalty <- -0.5*((kurt[-Var.xsim]-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
# penalty <- sum(penalty)

  # Simplistic error control;  -800 is effectively like -Inf:
  if(is.infinite(llr) | is.na(llr)){llr <- -800}

  # trustregion is the maximum value of llr that we actually trust.
  # So if llr>trustregion, return a value less than trustregion instead.
  #
  # Return the negative log-lik (for nloptr is a minimizer)
  # Actually return the positive log-lik for regular terms
  if (is.numeric(control.llik$MCMLE.trustregion) && llr>control.llik$MCMLE.trustregion) {
    return(2*control.llik$MCMLE.trustregion - llr)
  } else {
    return(llr)
  }
}
#' @export
llik.grad.Kpenalty <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){

  Kpenalty <- grep("Taper_Penalty",colnames(xsim), fixed=TRUE)

  # Obtain canonical parameters incl. offsets and difference from sampled-from
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  # Calculate log-importance-weights
  basepred <- xsim %*% etaparam + lrowweights(xsim)
  
  # Calculate the estimating function values sans offset
  llg <- - lweighted.mean(xsim, basepred)
  llg <- t(ergm.etagradmult(theta, llg, etamap))
  
  maxbase <- max(basepred)
  lwi <- basepred - (maxbase + log(sum(exp(basepred-maxbase))) )
#
  logm2 <- log(colMeans(xsim[,-Kpenalty]^2))
  logm4 <- log(colMeans(xsim[,-Kpenalty]^4))
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
  logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
  kurt <- exp(logm4-2*logm2)
  kurt[is.nan(kurt) | is.na(kurt)] <- 3

  Ek <- colSums(sweep(xsim,1,exp(lwi),"*"))
  t2 <- sweep(sweep(xsim,2,Ek,"-"),1,exp(lwi),"*")
  g_M2  <- t(xsim[,-Kpenalty]^2) %*% t2
  g_M4  <- t(xsim[,-Kpenalty]^4) %*% t2
  g <- kurt * (g_M4-2*g_M2)
browser()
  llg <- llg + g + control.llik$MCMLE.kurtosis.penalty
  llg[is.na(llg) | is.nan(llg)] <- 0

  # Return negative grad
  llg
}

#' @export
eq.jacobian.Kpenalty <- function(theta, xsim, xsim.obs=NULL,
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
llik.hessian.Kpenalty <- function(theta, xsim, xsim.obs=NULL,
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
llik.hessian.Kpenalty.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){
H <-  numDeriv::hessian(llik.fun.Kpenalty,theta,
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
H
}

#' @export
llik.hessian.Kpenalty.direct <- function(theta, xsim, xsim.obs=NULL,
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
#     message("Kpenalty=")
#     print(kurt[-exterms])
#   }
    penalty
  }

  if(F){
   hessian.penalty <- numDeriv::hessian(kfn, theta, 
    method.args=list(eps=1e-2, d=0.1, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE),
    xsim=xsim, etamap=etamap, eta0=eta0, control.llik=control.llik)
   #print(dim(H))
   #print(dim(hessian.penalty))
   print((H))
   print((hessian.penalty))
    H+hessian.penalty
  }else{
    H
  }
}

#llik.grad.Kpenalty.orig <- llik.grad.IS
#llik.hessian.Kpenalty.orig <- llik.hessian.IS
# Numerical versions
#' @export
llik.grad.Kpenalty.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){
  numDeriv::grad(llik.fun.Kpenalty,theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                 )
}
#' @export
eq.jacobian.Kpenalty.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=control.ergm.tapered.loglik()
                     ){
# numDeriv::jacobian(eq.fun.Kpenalty,theta,
  nloptr::nl.jacobian(x0=theta,fn=eq.fun.Kpenalty,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                 )
}
