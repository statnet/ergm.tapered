#' @import ergm statnet.common network
InitErgmTerm.Stereo <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "coef", "m"),
                      vartypes = c("formula", "numeric", "numeric"),
                      defaultvalues = list(NULL, NULL, NULL),
                      required = c(TRUE, TRUE, TRUE))
  beta <- a$coef
  nws <- a$m

  m <- ergm_model(a$formula, nw, response=response, ...)
  NVL(nws) <- summary(m, nw, response=response)

  if(!is.null(beta)){ stereo.mult <- beta }else{ stereo.mult <- 1 }
  # TODO: Names matching here?
  if(length(stereo.mult)!=1){
     stop("Invalid stereoing parameter vector beta: ",
          "wrong number of parameters: expected ",
          "1, but got ",length(stereo.mult),".")
  }

  # Should be empty network statistics
  gs0 <- summary(m, NULL, response=response)

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

  cnt <- c(ergm_mk_std_op_namewrap(paste0('Stereo(',beta,')'))(param_names(m, canonical=TRUE)), "Stereo_Penalty")

  list(name="stereo_term", coef.names = cnt,
       inputs=c(beta, nws), # Note: what gets passed is the difference between the empty network and the observed network.
       auxiliaries = ~.submodel_and_summary(a$formula),
       dependence=TRUE, emptynwstats = c(gs0, -2*log(beta[1]*beta[1]+sum((gs0-nws)^2))),
       map = map, gradient = gradient, params = params, minpar=m$etamap$mintheta, maxpar=m$etamap$maxtheta)
}
