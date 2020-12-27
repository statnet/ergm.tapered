#' @import ergm statnet.common network
InitErgmTerm.M4 <- function(nw, arglist, response=NULL, ...){
  a <- ergm::check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "coef", "m"),
                      vartypes = c("formula", "numeric", "numeric"),
                      defaultvalues = list(NULL, NULL, NULL),
                      required = c(TRUE, TRUE, TRUE))
  f <- a$formula
  beta <- a$coef
  nws <- a$m
  if(length(f)==2) f <- statnet.common::nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm::ergm_model(f, nw,...)
  NVL(nws) <- summary(m)

  if(!is.null(beta)){ taper.mult <- beta }else{ taper.mult <- 1 }
  # TODO: Names matching here?
  if(length(nws)==length(taper.mult) & length(nws) > 0) {
    #message("Using a tapered version of the model (based on passed tapering scale).")
    beta<-taper.mult
  }else{if(length(taper.mult)==1){
    #message("Using a tapered version of the model (based on default tapering scale).")
    beta<-taper.mult / ((2^2) * nws)
  }else{
     stop("Invalid M4 parameter vector beta: ",
          "wrong number of parameters: expected ",
          length(nws),
          " or 1, but got ",length(taper.mult),".")
  }}
  beta<-rep(1,length(beta))

  inputs <- ergm::to_ergm_Cdouble(m)
  
  # Should be empty network statistics
  gs0 <- summary(m)

  map <- function(x, n, ...){
#   c(ergm.eta(x, m$etamap))
#   c(-exp(x))
    x
  }

  gradient <- function(x, n, ...){
#   cbind(ergm.etagrad(x, m$etamap))
#   cbind(c(1,0,0,0),c(0,1,0,0),c(0,0,-exp(x[3]),0),c(0,0,0,-exp(x[4])))
#   cbind(c(-exp(x[1]),0),c(0,-exp(x[2])))
#   diag(-exp(x),ncol=length(x))
    diag(length(x))
  }

  cnt <- c(paste0( 'M4(',param_names(m, canonical=FALSE),')'))
  params <- rep(list(NULL), nparam(m))
  names(params) <- cnt

# cnt <- c(paste0('M4(',param_names(m, canonical=FALSE),",",beta,')'))
  list(name="m4_term", coef.names = cnt,
       inputs=c(beta, inputs, gs0-nws), # Note: what gets passed is the difference between the empty network and the observed network.
       dependence=TRUE, emptynwstats = ((gs0-nws)^4)*beta)
      #map = map, gradient = gradient, params = params, minpar=m$etamap$mintheta, maxpar=m$etamap$maxtheta)
}
