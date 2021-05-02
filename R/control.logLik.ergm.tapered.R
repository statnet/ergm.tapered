#  File R/control.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
#' Auxiliary for Controlling Tapered logLik.ergm 
#' 
#' Auxiliary function as user interface for fine-tuning logLik.ergm algorithm,
#' which approximates log likelihood values.
#' 
#' This function is only used within a call to the \code{\link{logLik.ergm}}
#' function.
#' 
#' @templateVar MCMCType MCMC
#'
#' @param MCMLE.kurtosis.prior Logical: If TRUE, use a prior for the kurtosis.
#' @param MCMLE.kurtosis.location numeric The mean of the prior for kurtosis and/or the 
#' taget kurtosis value.
#' @param MCMLE.kurtosis.scale numeric The standard deviation of the prior for kurtosis.
#' @param MCMLE.kurtosis.penalty numeric: The size of the penalty for larger values of the tapering parameter on the log-likelihood.
#' @param MCMLE.kurtosis.equality Logical: If TRUE, constrain the kurtosis in expectation.
#' @param MCMLE.metric Method to calculate the loglikelihood approximation.
#' See [control.ergm()] for the standard [ergm()] options. The tapering models add
#' \code{} and .
#' @param \dots Additional arguments, passed to other functions This argument
#' is helpful because it collects any control parameters that have been
#' deprecated; a warning message is printed in case of deprecated arguments.
#' @return A list with arguments as components. It includes the [control.ergm()] arguments.
#' @seealso [control.ergm()].
#' @references \itemize{ 
#' * Fellows, I. and M. S. Handcock (2017), 
#' Removing Phase Transitions from Gibbs Measures. Volume 54 of 
#' Proceedings of Machine Learning Research, Fort Lauderdale,
#' FL, USA, pp. 289–297. PMLR.
#' 
#' }
#' @keywords models
#' @export control.ergm.tapered
control.logLik.ergm.tapered<-function(
                       MCMLE.metric=c("lognormal", "logtaylor",
                         "Median.Likelihood",
                         "EF.Likelihood", "naive",
                         "Kpenalty"),
                       MCMC.esteq.exclude.statistics=NULL,
                       MCMLE.kurtosis.prior=TRUE,
                       MCMLE.kurtosis.location=3.0,
                       MCMLE.kurtosis.scale=0.3,
                       MCMLE.kurtosis.penalty=2.0,
                       MCMLE.kurtosis.equality=FALSE,
                       ...
                       ){
  formal.args<-formals(sys.function())
  control <- control.logLik.ergm(...)
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  match.arg.pars <- c("MCMLE.metric")
  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))

  set.control.class("control.logLik.ergm")
}
