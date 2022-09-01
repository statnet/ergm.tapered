################################################################################

#' @templateVar name Taper
#' @title Taper 
#' @description 
#'    Adds the terms specified in \code{formula} to the model \emph{and}
#'    imposes the quadratic penalty of Fellows and Handcock
#'    (2017). 
#' @usage
#' # binary: Taper(formula=NULL, coef=NULL, m=NULL)
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
#' @concept tapered
#' @import ergm statnet.common network stats
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
