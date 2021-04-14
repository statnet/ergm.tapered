#' @import ergm statnet.common network
InitErgmTerm.Kpenalty <- function(nw, arglist, response=NULL, ...){
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
    c(ergm.eta(x[-n], m$etamap), -3*exp(x[n])/(1+exp(x[n])))
  }

  gradient <- function(x, n, ...){
    a <- ergm.etagrad(x[-n], m$etamap)
    a <- rbind(a,0)
    cbind(a, rep(c(0,-3*exp(x[n])/((1+exp(x[n]))^2)),c(nrow(a)-1,1)))
  }

# mintheta <- c(m$etamap$mintheta,-Inf)
# maxtheta <- c(m$etamap$maxtheta,0)
  cnt <- c(param_names(m, canonical=FALSE), "Taper_Penalty")

  params <- c(rep(list(NULL), nparam(m)),-log(2))
# params <- c(as.list(beta),-log(2))
# params <- rep(list(NULL), nparam(m)+1)
  names(params) <- c(param_names(m, canonical=FALSE),"Taper_Penalty")

  list(name="Kpenalty_term", coef.names = cnt,
       inputs=c(beta, nws),
       auxiliaries = ~.submodel_and_summary(a$formula),
       dependence=TRUE, emptynwstats = c(gs0, sum((gs0-nws)^2*beta)),
       map = map, gradient = gradient,
#      minpar=mintheta, maxpar=maxtheta,
       params = params)
}
