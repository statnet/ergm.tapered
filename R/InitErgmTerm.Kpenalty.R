################################################################################

#' @templateVar name Kpenalty
#' @title Kpenalty 
#' @description 
#'    Adds \emph{only} the quadratic penalty of Fellows and Handcock
#'    (2017) for the terms specified in \code{formula} to the model.
#' @usage
#' # binary: Kpenalty(formula=NULL, coef=NULL, m=NULL)
#'
#' @param formula a valid formula for a  standard ERGM
#' @param coef a numeric vector of coefficients giving the penalty coefficients \eqn{\beta} for the tapering of the terms.
#'   If \code{NULL} is passed, the tapering coefficients are set to \code{1/(4*summary(formula))}, the default
#'   in Fellows and Handcock (2017).
#'   If a numeric vector is given, there are interpreted as the tapering coefficients of the terms in the
#'   model, including the terms enclosed in \code{offset()}.
#'   If a numeric scalar is given, it is interpreted as a multiplier of the default tapering coefficients , that is,
#' @param m numeric vector. If given, is
#'    the value of the network statistic relative to which the model is
#'    tapered. If omitted, it defaults to that of the model's LHS
#'    network if \code{formula} is one-sided and that of the network on
#'    the LHS of \code{formula} if it is two-sided.
#'
#' @template ergmTerm-general
#'
#' @concept operator
#' @concept tapered
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
    blim <- c(3,3,1) # max= blim[2], min = 3/(1+2^blim[2])
    b <- blim[3]/(x[n]*x[n])
    r <- blim[2]*exp(log(2)*(b-blim[1]))/(1+exp(log(2)*(b-blim[1])))
    if(is.na(r) | is.infinite(r) | is.nan(r)) r <- blim[2]
    c(ergm.eta(x[-n], m$etamap), r)
  }

  gradient <- function(x, n, ...){
    a <- ergm.etagrad(x[-n], m$etamap)
    a <- rbind(a,0)
    blim <- c(3,3,1) # 0.0001, max= blim[2], min = 3/(1+2^blim[2])
    b <- blim[3]/(x[n]*x[n])
    r <- -2*blim[3]*blim[2]*log(2)*exp(log(2)*(b-blim[1]))/(b^1.5*(1+exp(log(2)*(b-blim[1])))^2)
    if(is.na(r) | is.infinite(r) | is.nan(r)) r <- 0
    cbind(a, rep(c(0,r),c(nrow(a)-1,1)))
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
