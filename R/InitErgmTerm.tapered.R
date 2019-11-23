#' @import ergm statnet.common network
InitErgmTerm.Taper <- function(nw, arglist, response=NULL, ...){
  a <- ergd::check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "coef", "m"),
                      vartypes = c("formula", "numeric", "numeric"),
                      defaultvalues = list(NULL, NULL, NULL),
                      required = c(TRUE, TRUE, TRUE))
  f <- a$formula
  beta <- a$coef
  nws <- a$m
  if(length(f)==2) f <- statnet.common::nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergd::ergm_model(f, nw,...)
  NVL(nws) <- summary(m)

  if(!is.null(beta)){ taper.mult <- beta }else{ taper.mult <- 1 }
  # TODO: Names matching here?
  if(length(nws)==length(taper.mult) & length(nws) > 1) {
    #message("Using a tapered version of the model (based on passed tapering scale).")
    beta<-taper.mult
  }else{if(length(taper.mult)==1){
    #message("Using a tapered version of the model (based on default tapering scale).")
    beta<-taper.mult / ((2^2) * nws)
  }else{
     stop("Invalid tapering parameter vector beta: ",
          "wrong number of parameters: expected ",
          length(nws),
          " or 1, but got ",length(taper.mult),".")
  }}

  inputs <- ergd::to_ergm_Cdouble(m)
  
  # Should be empty network statistics
  gs0 <- summary(m)

  map <- function(x, n, ...){
    c(ergm.eta(x, m$etamap), -1)
  }

  gradient <- function(x, n, ...){
    cbind(ergm.etagrad(x, m$etamap), 0)
  }

  cnt <- c(paste0('Taper(',param_names(m, canonical=TRUE),",",beta,')'), "Taper_Penalty")

  params <- rep(list(NULL), nparam(m))
  names(params) <- param_names(m, canonical=FALSE)

  list(name="taper_term", coef.names = cnt,
       inputs=c(beta, inputs, gs0-nws), # Note: what gets passed is the difference between the empty network and the observed network.
       dependence=TRUE, emptynwstats = c(gs0, sum((gs0-nws)^2*beta)),
       map = map, gradient = gradient, params = params, minpar=m$etamap$mintheta, maxpar=m$etamap$maxtheta)
}
