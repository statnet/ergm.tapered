
InitErgmTerm.Taper <- function(nw, arglist, response=NULL, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "coef", "m"),
                      vartypes = c("formula", "numeric", "numeric"),
                      defaultvalues = list(NULL, NULL, NULL),
                      required = c(TRUE, TRUE, FALSE))
  f <- a$formula
  beta <- a$coef
  if(length(f)==2) f <- nonsimp.update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)

  m <- ergm.getmodel(f, nw, response=response,...)
  NVL(nws) <- ergm.getglobalstats(nw, m, response=response)
  beta <- rep(beta, length.out=length(nws))
  Clist <- ergm.Cprepare(nw, m, response=response)

  inputs <- pack.Clist_as_num(Clist)
  
  gs0 <- ergm.emptynwstats.model(m)

  map <- function(x, n, ...){
    c(ergm.eta(x, m$etamap), -1)
  }

  gradient <- function(x, n, ...){
    cbind(ergm.etagrad(x, m$etamap), 0)
  }

  params <- rep(list(NULL), coef.length.model(m))
  names(params) <- coef.names.model(m, canonical=FALSE)

  cnt <- c(paste0('Taper(',m$coef.names,",",beta,')'), "Taper_Penalty")
  
  list(name="taper_term", coef.names = cnt,
       inputs=c(beta, inputs, gs0-nws), # Note: what gets passed is the difference between the empty network and the observed network.
       dependence=TRUE, emptynwstats = c(gs, sum((gs0-nws)^2*beta)),
  map = map, gradient = gradient, params = params)
}
