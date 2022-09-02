\name{ergm-terms-taper}
\alias{ergm-terms-taper}
\alias{Taper}
\alias{Stereo}
\alias{Kpenalty}
\alias{Var}
\alias{M4}
\docType{methods}
\title{Terms Implemented in the \code{ergm.tapered} Package}
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
      in Fellows and Handcock (2017).
      If a numeric vector is given, there are interpreted as the tapering coefficients of the terms in the
      model, including the terms enclosed in \code{offset()}.
      If a numeric scalar is given, it is interpreted as a multiplier of the default tapering coefficients , that is, \code{1/(4*summary(formula))}
      \code{m}, if given, is
      the value of the network statistic relative to which the model is
      tapered. If omitted, it defaults to that of the model's LHS
      network if \code{formula} is one-sided and that of the network on
      the LHS of \code{formula} if it is two-sided.
      
    }
    \item{\code{Stereo(formula, coef, m)} (binary) (operator)}{\emph{Stero Tapered ERGM terms}:

      Adds the terms specified in \code{formula} to the model \emph{and}
      imposes the Stereo penalty of Blackburn (2021). This is akin to using an inverse 
      stereographic projection onto a sphere. See Section 5.3 for a development.
      The stereo penalty is weaker than the variance tapering of Fellows and Handcock
      (2017). 

      \code{coef} is a numeric vector of coefficients giving the penalty coefficients \eqn{\beta} for the tapering of the terms.
      If \code{NULL} is passed, the tapering coefficients are set to \code{1/(4*summary(formula))}, the default
      in Fellows and Handcock (2017).
      If a numeric vector is given, there are interpreted as the tapering coefficients of the terms in the
      model, including the terms enclosed in \code{offset()}.
      If a numeric scalar is given, it is interpreted as a multiplier of the default tapering coefficients , that is, \code{1/(4*summary(formula))}
      \code{m}, if given, is
      the value of the network statistic relative to which the model is
      tapered. If omitted, it defaults to that of the model's LHS
      network if \code{formula} is one-sided and that of the network on
      the LHS of \code{formula} if it is two-sided.
      
    }
    \item{\code{Kpenalty(formula, coef, m)} (binary) (operator)}{\emph{Stereo Tapered ERGM terms}:

       Adds \emph{only} the quadratic penalty of Fellows and Handcock
      (2017) for the terms specified in \code{formula} to the model.

      \code{coef} is a numeric vector of coefficients giving the penalty coefficients \eqn{\beta} for the tapering of the terms.
      If \code{NULL} is passed, the tapering coefficients are set to \code{1/(4*summary(formula))}, the default
      in Fellows and Handcock (2017).
      If a numeric vector is given, there are interpreted as the tapering coefficients of the terms in the
      model, including the terms enclosed in \code{offset()}.
      If a numeric scalar is given, it is interpreted as a multiplier of the default tapering coefficients , that is, \code{1/(4*summary(formula))}
      \code{m}, if given, is
      the value of the network statistic relative to which the model is
      tapered. If omitted, it defaults to that of the model's LHS
      network if \code{formula} is one-sided and that of the network on
      the LHS of \code{formula} if it is two-sided.
      
    }
    \item{\code{Var(formula, coef, m)} (binary) (operator)}{\emph{Var ERGM terms}:

      Adds the variance of the terms specified in \code{formula} to the model. 
      This is sometimes applied directly to control the variance of the standard
      statistics, but more typically as part of a Tapered ERGM.
      \code{coef} is a numeric vector of coefficients giving the penalty coefficients \eqn{\beta} for the tapering of the terms.
      If \code{NULL} is passed, the tapering coefficients are set to \code{1/(4*summary(formula))}, the default
      in Fellows and Handcock (2017).
      If a numeric vector is given, there are interpreted as the tapering coefficients of the terms in the
      model, including the terms enclosed in \code{offset()}.
      If a numeric scalar is given, it is interpreted as a multiplier of the default tapering coefficients , that is, \code{1/(4*summary(formula))}
      \code{m}, if given, is
      the value of the network statistic relative to which the model is
      tapered. If omitted, it defaults to that of the model's LHS
      network if \code{formula} is one-sided and that of the network on
      the LHS of \code{formula} if it is two-sided.
      
    }
    \item{\code{M4(formula, coef, m)} (binary) (operator)}{\emph{M4 Tapered ERGM terms}:

      Adds the fourth moment of each statistics specified by the terms in \code{formula} to the model. 
      This is sometimes applied directly to control the fourth moment of the standard
      statistics, but more typically as part of a Tapered ERGM.

      \code{coef} is a numeric vector of coefficients giving the penalty coefficients \eqn{\beta} for the tapering of the terms.
      If \code{NULL} is passed, the tapering coefficients are set to \code{1/(4*summary(formula))}, the default
      in Fellows and Handcock (2017).
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
    \item{Blackburn, B and Handcock, MS (2022). Practical Network Modeling via Tapered
      Exponential-family Random Graph Models. \emph{Journal of Computational and Graphical Statistics}
      \doi{10.1080/10618600.2022.2116444}.}
  }
}
\seealso{\code{\link[=ergm-package]{ergm}} package, \code{\link{search.ergmTerms}}, \code{\link{ergm}}, \code{\link{network}}, \code{\link{\%v\%}}, \code{\link{\%n\%}}}
\examples{
data(florentine)
m <- summary(flomarriage~edges + triangles)
summary(ergm(flomarriage ~ Taper(~edges + triangles, coef=c(.1, .1),m=m)))
}
\keyword{models}
\concept{tapered}
