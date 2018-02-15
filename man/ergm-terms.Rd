\name{ergm-terms}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
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
}

\section{Terms}{
  \describe{
    \item{\code{Taper(formula, coef, m)} (binary) (operator)}{\emph{Tapered ERGM terms}:

      Adds the terms specified in \code{formula} to the model \emph{and}
      imposes the quadratic penalty of Fellows and Handcock
      (2017). \code{coef} is a numeric vector (recycled as needed)
      giving the penalty coefficients \eqn{\beta}. \code{m}, if given, is
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
summary(ergm(flomarriage ~ Taper(~edges + triangles, c(.1, .1))))
}
\keyword{models}
