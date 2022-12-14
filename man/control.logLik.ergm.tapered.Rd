% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control.logLik.ergm.tapered.R
\name{control.logLik.ergm.tapered}
\alias{control.logLik.ergm.tapered}
\title{Auxiliary for Controlling Tapered logLik.ergm}
\usage{
control.logLik.ergm.tapered(
  MCMLE.metric = c("lognormal", "logtaylor", "Median.Likelihood", "EF.Likelihood",
    "naive", "Kpenalty"),
  MCMC.esteq.exclude.statistics = NULL,
  MCMLE.kurtosis.location = 3,
  MCMLE.kurtosis.scale = 0.3,
  MCMLE.kurtosis.penalty = 2,
  ...
)
}
\arguments{
\item{MCMLE.metric}{Method to calculate the loglikelihood approximation.}

\item{MCMC.esteq.exclude.statistics}{vector vector of names of statistics to exclude from the estimating equations.}

\item{MCMLE.kurtosis.location}{numeric The mean of the prior for kurtosis and/or the
taget kurtosis value.}

\item{MCMLE.kurtosis.scale}{numeric The standard deviation of the prior for kurtosis.}

\item{MCMLE.kurtosis.penalty}{numeric: The size of the penalty for larger values of the tapering parameter on the log-likelihood.
See \code{\link[=control.ergm]{control.ergm()}} for the standard \code{\link[=ergm]{ergm()}} options. The tapering models add
\code{} and .}

\item{\dots}{Additional arguments, passed to other functions This argument
is helpful because it collects any control parameters that have been
deprecated; a warning message is printed in case of deprecated arguments.}
}
\value{
A list with arguments as components. It includes the \code{\link[=control.ergm]{control.ergm()}} arguments.
}
\description{
Auxiliary function as user interface for fine-tuning logLik.ergm algorithm,
which approximates log likelihood values.
}
\details{
This function is only used within a call to the \code{\link{logLik.ergm}}
function.
}
\references{
\itemize{
\itemize{
\item Fellows, I. and M. S. Handcock (2017),
Removing Phase Transitions from Gibbs Measures. Volume 54 of
Proceedings of Machine Learning Research, Fort Lauderdale,
FL, USA, pp. 289–297. PMLR.
\item Blackburn, B. and M. S. Handcock (2022),
Practical Network Modeling via Tapered Exponential-family Random Graph Models.
Journal of Computational and Graphical Statistics
\doi{10.1080/10618600.2022.2116444}.
}

}
}
\seealso{
\code{\link[=control.ergm]{control.ergm()}}.
}
\keyword{models}
