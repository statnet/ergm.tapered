#' @export
llik.fun.obs.Kpenalty <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=ergm::control.ergm.tapered.loglik()
                     ){

  Kpenalty <- grep("Taper_Penalty",colnames(xsim), fixed=TRUE)

# Convert theta to eta
  eta <- ergm.eta(theta, etamap)
  etaparam <- eta-eta0

  basepred <- xsim %*% etaparam
  obspred <- xsim.obs %*% etaparam
  mb <- lweighted.mean(basepred,rowweights(xsim))
  vb <- lweighted.var(basepred,rowweights(xsim))
  mm <- lweighted.mean(obspred,lrowweights(xsim.obs))
  vm <- lweighted.var(obspred,lrowweights(xsim.obs))

# This is log-lik (and not its negative)
  llr <- (mm + control.llik$MCMLE.varweight*vm) - (mb + control.llik$MCMLE.varweight*vb)
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
  penalty <- -0.5*((kurt-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
  Kurt_penalty0 <- sum(penalty)

  maxbase <- max(basepred)
  lwi <- basepred - (maxbase + log(sum(exp(basepred-maxbase))) )
  logm2 <- log(colSums(sweep(xsim[,-Kpenalty]^2  ,1,exp(lwi),"*")))
  logm4 <- log(colSums(sweep(xsim[,-Kpenalty]^4,1,exp(lwi),"*")))
  logm2[is.infinite(logm2) | is.nan(logm2) | is.na(logm2)] <- 0
  logm4[is.infinite(logm4) | is.nan(logm4) | is.na(logm4)] <- 0
  kurt <- exp(logm4-2*logm2)
  kurt[is.nan(kurt) | is.na(kurt)] <- 3
  penalty <- -0.5*((kurt-control.llik$MCMLE.kurtosis.location)/control.llik$MCMLE.kurtosis.scale)^2
  Kurt_penalty <- sum(penalty)
#
  Kurt_penalty_diff <- Kurt_penalty - Kurt_penalty0
  Tpenalty_diff <- eta[Kpenalty] - eta0[Kpenalty]
#
  llr <- llr + 2*Tpenalty_diff + Kurt_penalty_diff

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

# Numerical versions
#' @export
llik.hessian.obs.Kpenalty.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=ergm::control.ergm.tapered.loglik()
                     ){
  numDeriv::hessian(llik.fun.obs.Kpenalty,theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                 )
}

#' @export
llik.grad.obs.Kpenalty.numDeriv <- function(theta, xsim, xsim.obs=NULL,
                     eta0, etamap, 
                     control.llik=ergm::control.ergm.tapered.loglik()
                     ){
  numDeriv::grad(llik.fun.obs.Kpenalty,theta,
                 xsim=xsim, xsim.obs=xsim.obs,
                 eta0=eta0, etamap=etamap,
                 control.llik=control.llik
                 )
}
