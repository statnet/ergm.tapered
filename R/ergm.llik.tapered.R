# Kurtosis related functions
eq.fun.tapered <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=list(MCMLE.trustregion=20, MCMLE.varweight=0.5, MCMLE.dampening=FALSE,MCMLE.dampening.min.ess=100, MCMLE.dampening.level=0.1, MCMLE.kurtosis.location=3.0, MCMLE.kurtosis.scale=0.3, MCMLE.kurtosis.penalty=2, MCMLE.kurtosis.prior=FALSE)
                     ){

  Ki <- control.llik$MCMLE.kurtosis.location

  eta <- ergm.eta(theta, etamap)
  eta[is.na(eta)] <- 0

  etaparam <- (eta-eta0)
  basepred <- xsim %*% etaparam
  maxbase <- max(basepred)
  lwi <- basepred - (maxbase + log(sum(exp(basepred-maxbase))) )
  Ek <- colSums(sweep(xsim,1,exp(lwi),"*"))
# logm2 <- log(colSums(sweep(xsim[,Var.xsim]  ,1,exp(lwi),"*")))
# logm4 <- log(colSums(sweep(xsim[,Var.xsim]^2,1,exp(lwi),"*")))
# logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
# logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
## kurt <- exp(logm4-2*logm2)
# kurt <- logm4-2*logm2
# kurt[is.nan(kurt) | is.na(kurt)] <- 3
# names(kurt) <- colnames(xsim)[Var.xsim]
# penalty <- -0.5*((kurt-control.llik$MCMLE.kurtosis.location)/0.6)^2
# penalty <- sum(penalty)
## cat(sprintf("penalty = %9.5g\n",penalty))
## print(round(kurt,2))
  Ek[-length(Ek)]
}
#end

llik.fun.tapered <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=list(MCMLE.trustregion=20, MCMLE.varweight=0.5, MCMLE.dampening=FALSE,MCMLE.dampening.min.ess=100, MCMLE.dampening.level=0.1, MCMLE.kurtosis.location=3.0, MCMLE.kurtosis.scale=0.3, MCMLE.kurtosis.penalty=2, MCMLE.kurtosis.prior=FALSE)
                     ){
# Convert theta to eta
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  basepred <- xsim %*% etaparam
  mb <- lweighted.mean(basepred,rowweights(xsim))
  vb <- lweighted.var(basepred,rowweights(xsim))
# This is log-lik (and not its negative)
  llr <- -mb - control.llik$MCMLE.varweight*vb

  # Simplistic error control;  -800 is effectively like -Inf:
  if(is.infinite(llr) | is.na(llr)){llr <- -800}

  # trustregion is the maximum value of llr that we actually trust.
  # So if llr>trustregion, return a value less than trustregion instead.
  #
  # Return the negative log-lik (for nloptr is a minimizer)
# if (is.numeric(control.llik$MCMLE.trustregion) && llr< -1000) {
#   return(1000+log(-llr))
# } else {
    return(-llr)
# }
}
llik.grad.tapered <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=list(MCMLE.trustregion=20, MCMLE.varweight=0.5, MCMLE.dampening=FALSE,MCMLE.dampening.min.ess=100, MCMLE.dampening.level=0.1, MCMLE.kurtosis.location=3.0, MCMLE.kurtosis.scale=0.3, MCMLE.kurtosis.penalty=2, MCMLE.kurtosis.prior=FALSE)
                     ){

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

eq.jacobian.tapered <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=list(MCMLE.trustregion=20, MCMLE.varweight=0.5, MCMLE.dampening=FALSE,MCMLE.dampening.min.ess=100, MCMLE.dampening.level=0.1, MCMLE.kurtosis.location=3.0, MCMLE.kurtosis.scale=0.3, MCMLE.kurtosis.penalty=2, MCMLE.kurtosis.prior=FALSE)
                     ){

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
# jac_Var <- t(xsim[,Var.xsim]) %*% t2
# jac_M4  <- t(xsim[,Var.xsim]^2) %*% t2
  jac   <-  t(xsim) %*% t2
# jac <- rbind(jac0, jac_M4-2*jac_Var)
  jac[is.na(jac) | is.nan(jac)] <- 0
  jac[-nrow(jac),]
}

llik.hessian.tapered <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=list(MCMLE.trustregion=20, MCMLE.varweight=0.5, MCMLE.dampening=FALSE,MCMLE.dampening.min.ess=100, MCMLE.dampening.level=0.1, MCMLE.kurtosis.location=3.0, MCMLE.kurtosis.scale=0.3, MCMLE.kurtosis.penalty=2, MCMLE.kurtosis.prior=FALSE)
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

llik.hessian.tapered.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=list(MCMLE.trustregion=20, MCMLE.varweight=0.5, MCMLE.dampening=FALSE,MCMLE.dampening.min.ess=100, MCMLE.dampening.level=0.1, MCMLE.kurtosis.location=3.0, MCMLE.kurtosis.scale=0.3, MCMLE.kurtosis.penalty=2, MCMLE.kurtosis.prior=FALSE)
                     ){
H <-  numDeriv::hessian(llik.fun.tapered,theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0, etamap, 
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
llik.grad.tapered.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=list(MCMLE.trustregion=20, MCMLE.varweight=0.5, MCMLE.dampening=FALSE,MCMLE.dampening.min.ess=100, MCMLE.dampening.level=0.1, MCMLE.kurtosis.location=3.0, MCMLE.kurtosis.scale=0.3, MCMLE.kurtosis.penalty=2, MCMLE.kurtosis.prior=FALSE)
                     ){
  numDeriv::grad(llik.fun.tapered,theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0, etamap, 
                 control.llik=control.llik
                 )
}
eq.jacobian.tapered.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=list(MCMLE.trustregion=20, MCMLE.varweight=0.5, MCMLE.dampening=FALSE,MCMLE.dampening.min.ess=100, MCMLE.dampening.level=0.1, MCMLE.kurtosis.location=3.0, MCMLE.kurtosis.scale=0.3, MCMLE.kurtosis.penalty=2, MCMLE.kurtosis.prior=FALSE)
                     ){
# numDeriv::jacobian(eq.fun.tapered,theta,
  nloptr::nl.jacobian(x0=theta,fn=eq.fun.tapered,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0, etamap, 
                 control.llik=control.llik
                 )
}
