#  File R/control.ergm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/citation
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
#' @param MCMLE.termination The criterion used for terminating MCMLE
#' estimation. The default for \code{ergm.tapered} is
#' `"precision"`: Terminate when the estimated loss in estimating precision
#' due to using MCMC standard errors is below the precision bound specified by
#' \code{MCMLE.MCMC.precision}, and the Hummel step length is 1 for two
#' consecutive iterations. See \code{MCMLE.MCMC.precision} for details. This
#' feature is in experimental status until we verify the coverage of the
#' standard errors. See the documentation for \code{\link{control.ergm}} for the other options
#' @param MCMLE.MCMC.precision
#' \code{MCMLE.MCMC.precision} is a vector of upper bounds on the standard
#' errors induced by the MCMC algorithm, expressed as a percentage of the total
#' standard error. The MCMLE algorithm will terminate when the MCMC standard
#' errors are below the precision bound, and the Hummel step length is 1 for
#' two consecutive iterations. This is an experimental feature [ergm()]. 
#' The default value in [ergm.tapered()] is 0.15, higher than in [ergm()].
#' @param estimate.tapered.bias logical. It \code{TRUE} the bias of the estimated estimates due to tapering is estimated.
#' If this is \code{FALSE} the MPLE of the untapered model is not computed, saving computational time. 
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
#' * Blackburn, B. and M. S. Handcock (2022), 
#' Practical Network Modeling via Tapered Exponential-family Random Graph Models.
#' Journal of Computational and Graphical Statistics
#' \doi{10.1080/10618600.2022.2116444}.
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
                       MCMLE.termination=c("precision","confidence", "Hummel", "Hotelling", "none"),
                       MCMLE.MCMC.precision=0.005,
                       estimate.tapered.bias=TRUE,
                 #     MCMC.effectiveSize.maxruns=8,
                       loglik=control.logLik.ergm.tapered(),
                       ...
                       ){
  formal.args<-formals(sys.function())
  control <- control.ergm(...)
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  match.arg.pars <- c("MCMLE.metric","MCMLE.termination")
  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))

  set.control.class(c("control.ergm", "control.ergm.tapered"))
}

STATIC_TAPERING_CONTROLS <- c("MCMLE.kurtosis.location", "MCMLE.kurtosis.scale", "MCMLE.kurtosis.penalty", "estimate.tapered.bias")
