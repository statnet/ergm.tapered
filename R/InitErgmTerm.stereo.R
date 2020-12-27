#' @import ergm statnet.common network
InitErgmTerm.Stereo <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "coef", "m"),
                      vartypes = c("formula", "numeric", "numeric"),
                      defaultvalues = list(NULL, NULL, NULL),
                      required = c(TRUE, TRUE, TRUE))
  f <- a$formula
  beta <- a$coef
  nws <- a$m
  if(length(f)==2) f <- statnet.common::nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm_model(f, nw, ...)
  NVL(nws) <- summary(m)

  if(!is.null(beta)){ stereo.mult <- beta }else{ stereo.mult <- 1 }
  # TODO: Names matching here?
  if(length(stereo.mult)!=1){
     stop("Invalid stereoing parameter vector beta: ",
          "wrong number of parameters: expected ",
          "1, but got ",length(stereo.mult),".")
  }

  inputs <- to_ergm_Cdouble(m)
  
  # Should be empty network statistics
  gs0 <- summary(m)

# This is a fudge line
  beta <-rep(beta, length(gs0))

  map <- function(x, n, ...){
    c(ergm.eta(x, m$etamap), -1)
  }

  gradient <- function(x, n, ...){
    cbind(ergm.etagrad(x, m$etamap), 0)
  }

  params <- rep(list(NULL), nparam(m))
  names(params) <- param_names(m, canonical=FALSE)

  cnt <- c(paste0('Stereo(',param_names(m, canonical=TRUE),",",beta,')'), "Stereo_Penalty")
  #print(beta)
  #print(nws)
  list(name="stereo_term", coef.names = cnt,
       inputs=c(beta, inputs, gs0-nws), # Note: what gets passed is the difference between the empty network and the observed network.
       dependence=TRUE, emptynwstats = c(gs0, -2*log(beta[1]*beta[1]+sum((gs0-nws)^2))),
       map = map, gradient = gradient, params = params, minpar=m$etamap$mintheta, maxpar=m$etamap$maxtheta)
}
