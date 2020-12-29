#' @import ergm statnet.common network
InitErgmTerm.Taper <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
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
  
  # Should be empty network statistics
  gs0 <- summary(m, NULL, response=response)

  map <- function(x, n, ...){
    c(ergm.eta(x, m$etamap), -1)
  }

  gradient <- function(x, n, ...){
    cbind(ergm.etagrad(x, m$etamap), 0)
  }

  cnt <- c(ergm_mk_std_op_namewrap(paste0('Taper(',beta,')'))(param_names(m, canonical=TRUE)), "Taper_Penalty")

  params <- rep(list(NULL), nparam(m))
  names(params) <- param_names(m, canonical=FALSE)

  list(name="taper_term", coef.names = cnt,
       inputs=c(beta, nws),
       auxiliaries = ~.submodel_and_summary(a$formula),
       dependence=TRUE, emptynwstats = c(gs0, sum((gs0-nws)^2*beta)),
       map = map, gradient = gradient, params = params, minpar=m$etamap$mintheta, maxpar=m$etamap$maxtheta)
}
