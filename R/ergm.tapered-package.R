#' Fit, Simulate and Diagnose Tapered Exponential-Family Models for Networks
#' 
#' \code{\link[=ergm.tapered-package]{ergm.tapered}} is a collection of functions to plot, fit,
#' diagnose, and simulate from Tapered exponential-family random graph models (ERGMs).
#' For a list of functions type: \code{help(package='ergm.tapered')}
#' 
#' A good place to start is the vignette at \url{https://github.com/statnet/ergm.tapered}
#' and the first two referenced papers below.
#
#' For a complete list of the functions, use \code{library(help="ergm.tapered")} or
#' read the rest of the manual. 
#' 
#' When publishing results obtained using this package, please cite the
#' original authors as described in \code{citation(package="ergm.tapered")}.
#' 
#' All programs derived from this package must cite it.
#' 
#' This package is founded on \code{\link[=ergm-package]{ergm}} and can be thought of as a generalization of it.
#' Most of its functionality is from \code{\link[=ergm-package]{ergm}}.  
#' 
#' The
#' \code{\link[=ergm-package]{ergm}} package implements maximum likelihood
#' estimates of ERGMs to be calculated using Markov Chain Monte Carlo (via
#' \code{\link{ergm}}). The package also provides tools for simulating networks
#' (via \code{\link{simulate.ergm}}) and assessing model goodness-of-fit (see
#' \code{\link{mcmc.diagnostics}} and \code{\link{gof.ergm.tapered}}).
#' 
#' For detailed information on how to download and install the software, go to
#' the \code{\link[=ergm.tapered-package]{ergm.tapered}} website: \url{https://statnet.org}. A
#' tutorial, support newsgroup, references and links to further resources are
#' provided there.
#' 
#' @name ergm.tapered-package
#' @docType package
#' @author Mark S. Handcock \email{handcock@stat.ucla.edu},\cr Pavel N. Krivitsky
#' \email{pavel@statnet.org}, and\cr Ian E. Fellows
#' \email{ian@fellstat.com}
#' 
#' Maintainer: Mark S. Handcock \email{handcock@@stat.ucla.edu}
#' @references \itemize{ 
#' * Fellows, I. and M. S. Handcock (2017), 
#' Removing Phase Transitions from Gibbs Measures. Volume 54 of 
#' Proceedings of Machine Learning Research, Fort Lauderdale,
#' FL, USA, pp. 289–297. PMLR.
#' * Blackburn, B. and M. S. Handcock (2022), 
#' Practical Network Modeling via Tapered Exponential-family Random Graph Models.
#' Journal of Computational and Graphical Statistics
#' \doi{10.1080/10618600.2022.2116444}.
#' * Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003a).
#' \pkg{ergm}: A Package to Fit, Simulate and Diagnose Exponential-Family
#' Models for Networks.  Statnet Project, Seattle, WA.  Version 3,
#' \url{https://statnet.org}.
#' * Handcock MS, Hunter DR, Butts CT, Goodreau SM, Morris M (2003b).
#' \pkg{statnet}: Software Tools for the Statistical Modeling of Network Data.
#' Statnet Project, Seattle, WA.  Version 3, \url{https://statnet.org}.
#' 
#' }
#' @keywords package models
#' @useDynLib ergm.tapered
NULL
