#  File R/control.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################
#' Auxiliary for Controlling Tapered ERGM Fitting
#' 
#' Auxiliary function as user interface for fine-tuning 'ergm' fitting
#' for tapered models.
#' 
#' This function is only used within a call to the [ergm.tapered()] function.
#' See the \code{usage} section in [ergm()] for details.
#' 
#' @templateVar MCMCType MCMC
#'
#' @param MCMLE.metric Method to calculate the loglikelihood approximation.
#' @param MCMC.esteq.exclude.statistics vector vector of names of statistics to exclude from the estimating equations.
#' @param MCMLE.kurtosis.location numeric The mean of the prior for kurtosis and/or the 
#' target kurtosis value.
#' @param MCMLE.kurtosis.scale numeric The standard deviation of the prior for kurtosis.
#' @param MCMLE.kurtosis.penalty numeric: The size of the penalty for larger values of the tapering parameter on the log-likelihood.
#' See [control.ergm()] for the standard [ergm()] options. The tapering models add
#' \code{} and .
#' @param loglik list List of additional control arguments for the loglik. See \code{\link{control.logLik.ergm.tapered}}.
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
control.ergm.tapered<-function(
                       MCMLE.metric=c("lognormal", "logtaylor",
                         "Median.Likelihood",
                         "EF.Likelihood", "naive",
                         "Kpenalty"),
                       MCMC.esteq.exclude.statistics=NULL,
                       MCMLE.kurtosis.location=3.0,
                       MCMLE.kurtosis.scale=0.3,
                       MCMLE.kurtosis.penalty=2.0,
                       loglik=control.logLik.ergm.tapered(),
                       ...
                       ){
  formal.args<-formals(sys.function())
  control <- control.ergm(...)
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  match.arg.pars <- c("MCMLE.metric")
  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))

  set.control.class(c("control.ergm", "control.ergm.tapered"))
}

STATIC_TAPERING_CONTROLS <- c("MCMLE.kurtosis.location", "MCMLE.kurtosis.scale", "MCMLE.kurtosis.penalty")
