# Version working on the tapering model with just a penalty
# for Kurtosis and large tau (or small r > 0)
#' @keywords internal 
#' @export
llik.fun.Kpenalty <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=ergm::control.ergm.tapered.loglik()
                     ){

  Kpenalty <- grep("Taper_Penalty",colnames(xsim), fixed=TRUE)
  blim <- c(3,3) # max= blim[2], min = 3/(1+2^blim[2])

# Convert theta to eta
  eta <- ergm.eta(theta, etamap)
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

  Tpenalty_diff <- eta[Kpenalty] - eta0[Kpenalty]

  if(eta[Kpenalty] < (blim[2]-0.001) | eta0[Kpenalty] < (blim[2]-0.001)){
  logm2 <- log(colMeans(xsim[,-Kpenalty]^2))
  logm4 <- log(colMeans(xsim[,-Kpenalty]^4))
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
  logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
  kurt <- exp(logm4-2*logm2)
  kurt[is.nan(kurt) | is.na(kurt)] <- blim[2]
  penalty <- -0.5*((kurt-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
  Kurt_penalty0 <- mean(penalty)

  maxbase <- max(basepred)
  lwi <- basepred - (maxbase + log(sum(exp(basepred-maxbase))) )
  logm2 <- log(colSums(sweep(xsim[,-Kpenalty]^2,1,exp(lwi),"*")))
  logm4 <- log(colSums(sweep(xsim[,-Kpenalty]^4,1,exp(lwi),"*")))
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
  logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
  kurt <- exp(logm4-2*logm2)

  kurt[is.nan(kurt) | is.na(kurt)] <- blim[2]
  penalty <- -0.5*((kurt-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
  Kurt_penalty <- mean(penalty)
#
  Kurt_penalty_diff <- Kurt_penalty - Kurt_penalty0
  llr <- llr + Kurt_penalty_diff
  }
  llr <- llr + control.llik$MCMLE.kurtosis.penalty*Tpenalty_diff

  # Simplistic error control;  -800 is effectively like -Inf:
  if(is.infinite(llr) | is.na(llr)){llr <- -800}

  # Return the negative log-lik (for nloptr is a minimizer)
  # Actually return the positive log-lik for regular terms
  llr
}

# Numerical versions
#' @keywords internal 
#' @export
llik.hessian.Kpenalty.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=ergm::control.ergm.tapered.loglik()
                     ){
  numDeriv::hessian(llik.fun.Kpenalty,theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                 )
}

#' @keywords internal 
#' @export
llik.grad.Kpenalty.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=ergm::control.ergm.tapered.loglik()
                     ){
  numDeriv::grad(llik.fun.Kpenalty,theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                )
}
#' @keywords internal 
#' @export
llik.grad.Kpenalty <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=ergm::control.ergm.tapered.loglik()
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

  logm2 <- log(colSums(sweep(xsim[,-Kpenalty]^2,1,exp(lwi),"*")))
  logm4 <- log(colSums(sweep(xsim[,-Kpenalty]^4,1,exp(lwi),"*")))
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
  logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
  kurt <- exp(logm4-2*logm2)
  kurt[is.nan(kurt) | is.na(kurt)] <- 3

  Ek <- colSums(sweep(xsim,1,exp(lwi),"*"))
  dwdeta <- sweep(sweep(xsim,2,Ek,"-"),1,exp(lwi),"*")
  g_M2  <- t(xsim[,-Kpenalty]^2) %*% dwdeta
  g_M4  <- t(xsim[,-Kpenalty]^4) %*% dwdeta
  dkdeta <- (g_M4*exp(2*logm2)-2*exp(logm4+logm2)*g_M2) / exp(4*logm2)
  g <- -((kurt-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)*dkdeta / control.llik$MCMLE.kurtosis.scale
  g <- colMeans(g)
#
   llik.fun.Kpenalty.tau <- function(tau, Kpenalty, theta ,xsim, xsim.obs,
                 eta0, etamap,
                 control.llik){
    theta[Kpenalty] <- tau
    llik.fun.Kpenalty(theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                )
   }
   aa=numDeriv::grad(llik.fun.Kpenalty.tau,theta[Kpenalty],Kpenalty=Kpenalty,theta=theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                )
  llg <- as.numeric(llg) + g
  llg[Kpenalty] <- aa
  llg[is.na(llg) | is.nan(llg)] <- 0

  # Return negative grad
  llg
}

#' @keywords internal 
#' @export
llik.hessian.Kpenalty <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=ergm::control.ergm.tapered.loglik()
                     ){

  Kpenalty <- grep("Taper_Penalty",colnames(xsim), fixed=TRUE)

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

  maxbase <- max(basepred)
  lwi <- basepred - (maxbase + log(sum(exp(basepred-maxbase))) )

  m2 <- colSums(sweep(xsim[,-Kpenalty]^2,1,exp(lwi),"*"))
  m4 <- colSums(sweep(xsim[,-Kpenalty]^4,1,exp(lwi),"*"))
  logm2 <- log(m2)
  logm4 <- log(m4)
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
  logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
  kurt <- exp(logm4-2*logm2)
  kurt[is.nan(kurt) | is.na(kurt)] <- 3

  Ek <- colSums(sweep(xsim,1,exp(lwi),"*"))
  xE <- sweep(xsim,2,Ek,"-")
  dwdeta <- sweep(xE,1,exp(lwi),"*")
  g_M2  <- t(xsim[,-Kpenalty]^2) %*% dwdeta
  g_M4  <- t(xsim[,-Kpenalty]^4) %*% dwdeta
  dkdeta <- (g_M4*exp(2*logm2)-2*exp(logm4+logm2)*g_M2) / exp(4*logm2)

  dgrad <- matrix(0,ncol=ncol(H),nrow=nrow(H))
  for(l in 1:length(kurt)){
   dgradkurt <- matrix(0,ncol=ncol(H),nrow=nrow(H))
   ddM2 <- t(sweep(xE,1,xsim[,l]^2-m2[l],"*")) %*% dwdeta
   ddM4 <- t(sweep(xE,1,xsim[,l]^4-m4[l],"*")) %*% dwdeta
   for(i in 1:nrow(dgradkurt)){
    for(j in 1:ncol(dgradkurt)){
     dKdj_n <- g_M4[l,j] * m2[l]*m2[l] - 2*m2[l]*m4[l]*g_M2[l,j] 
     dKddij <- ddM4[i,j] * m2[l]*m2[l]+g_M4[l,j]*2*m2[l]*g_M2[l,i] - 2*g_M2[l,j]*g_M2[l,i]*m4[l] - 2*m2[l]*(ddM2[i,j]*m4[l] + g_M4[l,i]*g_M2[l,j]) 
     dKddij <- dKddij * exp(4*logm2[l]) - dKdj_n*4*exp(3*logm2[l])*g_M2[l,i]
     dKddij <- dKddij / exp(8*logm2[l])
     dgradkurt[i,j] <- - dkdeta[l,i]*dkdeta[l,j]/(control.llik$MCMLE.kurtosis.scale*control.llik$MCMLE.kurtosis.scale) -
     ((kurt[l]-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)* dKddij / control.llik$MCMLE.kurtosis.scale
    }
   }
   dgradkurt[is.na(dgradkurt) | is.nan(dgradkurt) | is.infinite(dgradkurt)] <- 0
   dgrad <- dgrad + dgradkurt / length(kurt)
  }
  
   llik.grad.Kpenalty.tau <- function(tau, Kpenalty, theta ,xsim, xsim.obs,
                 eta0, etamap,
                 control.llik){
    theta[Kpenalty] <- tau
    llik.grad.Kpenalty(theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                )
   }
   aa=as.numeric(numDeriv::jacobian(llik.grad.Kpenalty.tau,theta[Kpenalty],Kpenalty=Kpenalty,theta=theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                ))
  HA <- H + dgrad
  HA[ Kpenalty,] <- aa
  HA[,Kpenalty ] <- aa

  HA[is.na(HA) | is.nan(HA) | is.infinite(HA)] <- 0
  # a simple check for p.s.d. of the approximation
  if(any(diag(HA) > 0)){
   H[is.na(H) | is.nan(H) | is.infinite(H)] <- 0
   H
  }else{
   HA
  }
}
