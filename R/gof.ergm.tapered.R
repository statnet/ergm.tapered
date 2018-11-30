#' Fits a Tapered ERGM
#' @param formula An ergm formula to fit
#' @param r The scaling factor to use for the hueristic of setting beta equal to r the standard deviation of the observed statistics
#' @param beta The tapering parameters. If not null, these override the hueristics.
#' @param tapering.centers The centers of the tapering terms. If null, these are taken to be the mean value parameters.
#' @param control An object of class control.ergm. Passed to the ergm function.
#' @param ... Additional arguments to ergm.
#' @returns
#' An object of class c('tapered.ergm','ergm') containing the fit model. In addition to all of the ergm items, 
#' this object contains tapering.centers, tapering.coef and orig.formula. tapering.centers are the centers for the tapering term.
#' tapering.coef are the tapering coefficients = 1/ beta^2. orig.formula is the formula passed into ergm.tapered.
#' @importFrom stats var as.formula
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
