################################################################################

#' @templateVar name Var
#' @title Var 
#' @description 
#'    Adds the variance of the terms specified in \code{formula} to the model. 
#'    This is sometimes applied directly to control the variance of the standard
#'    statistics, but more typically as part of a Tapered ERGM.
#' @usage
#' # binary: Var(formula=NULL, coef=NULL, m=NULL)
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
#' @template ergmTerm-taper-references
#'
#' @concept operator
#' @concept tapered
#' @import ergm statnet.common network
InitErgmTerm.Var <- function(nw, arglist, response=NULL, ...){
  a <- ergm::check.ErgmTerm(nw, arglist,
                      varnames = c("formula", "coef", "m"),
                      vartypes = c("formula", "numeric", "numeric"),
                      defaultvalues = list(NULL, NULL, NULL),
                      required = c(TRUE, FALSE, FALSE))
  beta <- a$coef
  nws <- a$m

  m <- ergm_model(a$formula, nw,...)
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
  beta<-rep(1,length(beta))

  # Should be empty network statistics
  gs0 <- summary(m, NULL, response=response)

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

  cnt <- c(paste0('Var(',param_names(m, canonical=FALSE),')'))
  params <- rep(list(NULL), nparam(m))
  names(params) <- cnt

# cnt <- c(paste0('Var(',param_names(m, canonical=FALSE),",",beta,')'))
  list(name="var_term", coef.names = cnt,
       inputs=c(nws),
       auxiliaries = ~.submodel_and_summary(a$formula),
       dependence=TRUE, emptynwstats = c((gs0-nws)^2*beta))
# The next is the key for curved version of this term
#      map = map, gradient = gradient, params = params, minpar=m$etamap$mintheta, maxpar=m$etamap$maxtheta)
}
