#' Conduct Goodness-of-Fit Diagnostics on a Tapered Exponential Family Random Graph
#' Model
#' 
#' \code{\link{gof}} calculates \eqn{p}-values for geodesic distance, degree,
#' and reachability summaries to diagnose the goodness-of-fit of tapered exponential
#' family random graph models.  See \code{\link{ergm.tapered}} for more information on
#' these models.
#' 
#' A sample of graphs is randomly drawn from the specified model.  The first
#' argument is typically the output of a call to \code{\link{ergm}} and the
#' model used for that call is the one fit.
#' 
#' For \code{GOF = ~model}, the model's observed sufficient statistics are
#' plotted as quantiles of the simulated sample. In a good fit, the observed
#' statistics should be near the sample median (0.5).
#'
#' @param object Either a formula or an \code{tapered.ergm} object.
#' See documentation for \code{\link{ergm.tapered}}.
#' An object of class c('tapered.ergm','ergm') contains the fit model. In addition to all of the ergm items, 
#' this object contains tapering.centers, tapering.coef and orig.formula. tapering.centers are the centers for the tapering term.
#' tapering.coef are the tapering coefficients = 1/ beta^2. orig.formula is the formula passed into ergm.tapered.
#' @param \dots Additional arguments, to be passed to lower-level functions.
#' These are the same as for \code{\link{gof.ergm}} and that documentation should be consulted.
#' @return \code{\link{gof}}, \code{\link{gof.ergm}}, and
#' \code{\link{gof.tapered.ergm}} return an object of class \code{tapered.ergm.gof}, which inherits from class `gof`.  This
#' is a list of the tables of statistics and \eqn{p}-values.  This is typically
#' plotted using \code{\link{plot.gof}}.
#' @seealso [ergm()], [network()], [simulate.ergm()], [summary.ergm()], [gof.ergm()]
#' @keywords models
#' @references 
#' Fellows, I and Handcock, MS (2017). Removing Phase Transitions from Gibbs Measures. Proceedings of Machine Learning Research, 54:289-297.
#' @examples 
#' \dontrun{
#' data(sampson)
#' fit <- ergm.tapered(samplike ~ edges + triangles())
#' summary(fit)
#' }
#' @export
gof.tapered.ergm <- function(object, ...){
  
  # do some formula magic
# .taper.coef <- object$tapering.coef
# .taper.center <- object$tapering.centers
# env <- new.env(parent=environment(object))
  env <- new.env(parent=emptyenv())
  env$.taper.center <- object$tapering.centers
  env$.taper.coef <- object$tapering.coef
  environment(object) <- env
  
  # GOF ergm
  class(object) <- "ergm"
  gof <- gof(object, ...)
  
  class(gof) <- c("tapered.ergm.gof", class(gof))
  
  gof
}
