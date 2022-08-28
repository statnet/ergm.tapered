\name{ergm-terms-taper}
\alias{ergm-terms-taper}
\alias{Taper}
\alias{InitErgmTerm.Taper}
\docType{methods}
\title{Terms Implemented in the \code{ergm.taper} Package}
\description{
  This page describes the possible terms (and hence network statistics)
  included in the \code{ergm.tapered} package.

  See the \code{\link[ergm]{ergm-terms}} documentation in the
  \code{ergm} package for a complete description of what ERGM terms are
  and how they are used.

  These terms are usually accessed via the \code{ergm.tapered} function call to mimic \code{ergm}.
  
}

\section{Terms}{
  \describe{
    \item{\code{Taper(formula, coef, m)} (binary) (operator)}{\emph{Tapered ERGM terms}:

      Adds the terms specified in \code{formula} to the model \emph{and}
      imposes the quadratic penalty of Fellows and Handcock
      (2017). \code{coef} is a numeric vector of coefficients giving the penalty coefficients \eqn{\beta} for the tapering of the terms.
      If \code{NULL} is passed, the tapering coefficients are set to \code{1/(4*summary(formula))}, the default
      in Fellows and Handcock (2016).
      If a numeric vector is given, there are interpreted as the tapering coefficients of the terms in the
      model, including the terms enclosed in \code{offset()}.
      If a numeric scalar is given, it is interpreted as a multiplier of the default tapering coefficients , that is, \code{1/(4*summary(formula))}
      \code{m}, if given, is
      the value of the network statistic relative to which the model is
      tapered. If omitted, it defaults to that of the model's LHS
      network if \code{formula} is one-sided and that of the network on
      the LHS of \code{formula} if it is two-sided.
      
    }
  }
}
\references{
  \itemize{
    \item{Fellows, I and Handcock, MS (2017). Removing Phase Transitions from
      Gibbs Measures. \emph{Proceedings of Machine Learning Research}, 54:289-297.}
  }
}
\seealso{\code{\link[=ergm-package]{ergm}} package, \code{\link{search.ergmTerms}}, \code{\link{ergm}}, \code{\link{network}}, \code{\link{\%v\%}}, \code{\link{\%n\%}}}
\examples{
data(florentine)
m <- summary(flomarriage~edges + triangles)
summary(ergm(flomarriage ~ Taper(~edges + triangles, coef=c(.1, .1),m=m)))
}
\keyword{models}
