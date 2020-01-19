#' @import ergd statnet.common network
InitErgmTerm.Kurt <- function(nw, arglist, response=NULL, ...){
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
  if(length(nws)==length(taper.mult) & length(nws) > 0) {
    #message("Using a tapered version of the model (based on passed tapering scale).")
    beta<-taper.mult
  }else{if(length(taper.mult)==1){
    #message("Using a tapered version of the model (based on default tapering scale).")
    beta<-taper.mult / ((2^2) * nws)
  }else{
     stop("Invalid Kurt parameter vector beta: ",
          "wrong number of parameters: expected ",
          length(nws),
          " or 1, but got ",length(taper.mult),".")
  }}
  beta<-rep(1,length(beta))

  inputs <- ergd::to_ergm_Cdouble(m)
  
  # Should be empty network statistics
  gs0 <- summary(m)

  map <- function(x, n, ...){
#   c(ergm.eta(x, m$etamap))
    c(6*exp(x),-exp(x)) 
  }

  gradient <- function(x, n, ...){
#   cbind(ergm.etagrad(x, m$etamap))
#  cbind(-6*diag(length(x)),diag(length(x)))
   cbind(6*diag(exp(x),ncol=length(x)),diag(-exp(x),ncol=length(x)))
  }

  cnt <- c(paste0('Var(',param_names(m, canonical=FALSE),')'),
           paste0( 'M4(',param_names(m, canonical=FALSE),')'))
  cnt_curved <- c(paste0('Kurt(',param_names(m, canonical=FALSE),')'))
  params <- rep(list(NULL), nparam(m))
  names(params) <- cnt_curved

# cnt <- c(paste0('Kurt(',param_names(m, canonical=FALSE),",",beta,')'))
  list(name="kurt_term", coef.names = cnt,
       inputs=c(beta, beta, inputs, gs0-nws), # Note: what gets passed is the difference between the empty network and the observed network.
       dependence=TRUE,  emptynwstats = c((gs0-nws)^2,(gs0-nws)^4),
       map = map, gradient = gradient, params = params, minpar=m$etamap$mintheta, maxpar=m$etamap$maxtheta)
}
