#' @import ergm statnet.common network
InitErgmTerm.Kurt <- function(nw, arglist, response=NULL, ...){
  a <- ergm::check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "coef", "m"),
                      vartypes = c("formula", "numeric", "numeric"),
                      defaultvalues = list(NULL, NULL, NULL),
                      required = c(TRUE, TRUE, TRUE))
  beta <- a$coef
  nws <- a$m

  m <- ergm_model(a$formula, nw, response=response, ...)
  NVL(nws) <- summary(m, nw, response=response)

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

  # Should be empty network statistics
  gs0 <- summary(m, NULL, response=response)

  map <- function(x, n, ...){
#   c(ergm.eta(x, m$etamap))
#   c(exp(x), -exp(3*x)*8/3)
    c((x), rep(0,length(x)))
  }

  gradient <- function(x, n, ...){
## cbind(-6*diag(length(x)),diag(length(x)))
## diag(-exp(x),ncol=length(x))
#  cbind(ergm.etagrad(x, m$etamap))
#  diag(length(x))
#  cbind(diag(exp(x),ncol=length(x)),diag(-exp(3*x)*8/2,ncol=length(x)))
   cbind(diag((x),ncol=length(x)),diag(rep(0,length(x)),ncol=length(x)))
  }

  cnt_Var <- paste0('Var(',param_names(m, canonical=FALSE),')')
  cnt <- c(paste0('Var(',param_names(m, canonical=FALSE),')'),
           paste0( 'M4(',param_names(m, canonical=FALSE),')'))
  cnt_curved <- c(paste0('Kurt(',param_names(m, canonical=FALSE),')'))
  params <- rep(list(NULL), nparam(m))
  names(params) <- cnt_Var

  beta <- rep(1,length(beta))
# cnt <- c(paste0('Kurt(',param_names(m, canonical=FALSE),",",beta,')'))
  list(name="kurt_term", coef.names = cnt,
#      inputs=c(beta, -8*beta^3/3, inputs, gs0-nws), # Note: what gets passed is the difference between the empty network and the observed network.
       inputs=c(beta, beta, nws), # Note: what gets passed is the difference between the empty network and the observed network.
       auxiliaries = ~.submodel_and_summary(a$formula),
       dependence=TRUE,  emptynwstats = c((gs0-nws)^2,(gs0-nws)^4),
       map = map, gradient = gradient, params = params, 
       minpar=m$etamap$mintheta,
       maxpar=m$etamap$maxtheta
               )
}
