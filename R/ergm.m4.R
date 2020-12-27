#' Fits a Tapered ERGM using kurtosis penalized likelihood
#' @param formula An ergm formula to fit
#' @param taper.terms Specification of the tapering used. If the character variable "dependence" then all the dependent
#' terms are tapered. If the character variable "all" then all terms are tapered.
#' It can also be the RHS of a formula giving the terms to be tapered. 
#' @param family The type of tapering used. This should either be the \code{stereo} or \code{taper}, the 
#' tapering model of Fellows and Handcock (2016).
#' @param r The (optional) scaling factor to use for the hueristic of setting the initial beta; equal to r standard deviations of the observed statistics. The default is 2.
#' @param beta The initial tapering parameters, expressed as in Fellows and Handcock (2017). If not NULL, these override the hueristics (r).
#' @param tau The initial tapering parameters, expressed as natural parameters. If not NULL, these override the beta and the hueristics (r).
#' @param tapering.centers The centers of the tapering terms. If null, these are taken to be the mean value parameters.
#' @param target.stats {vector of "observed network statistics,"
#' if these statistics are for some reason different than the 
#' actual statistics of the network on the left-hand side of
#' \code{formula}.
#' Equivalently, this vector is the mean-value parameter values for the
#' model.  If this is given, the algorithm finds the natural
#' parameter values corresponding to these mean-value parameters.
#' If \code{NULL}, the mean-value parameters used are the observed
#' statistics of the network in the formula.
#' }
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
#' fit <- ergm.kurtosis(samplike ~ edges + triangles(), tau=0.1))
#' summary(fit)
#' }
#' @export
ergm.m4 <- function(formula,
			 taper.terms="dependent",
			 family="taper",
                         r=2, beta=NULL, tau=NULL, tapering.centers=NULL, target.stats=NULL,
                         control = control.ergm(MCMLE.termination="confidence"), ...){
  
  # Determine the dyadic independence terms
  nw <- ergm.getnetwork(formula)
  m<-ergm_model(formula, nw)

  if(is.null(tapering.centers)) tapering.centers <- target.stats

  if(is.null(tapering.centers))
    ostats <- summary(formula)
  else
    ostats <- tapering.centers
  
  # set tapering terms
  if(is.character(taper.terms) & length(taper.terms)==1){
   if(taper.terms=="dependent"){
     a <- sapply(m$terms, function(term){is.null(term$dependence) || term$dependence})
     taper.terms <- list_rhs.formula(formula)
     tmp <- taper.terms
     taper.terms <- NULL
     for(i in seq_along(tmp)){if(a[i]){taper.terms <- c(taper.terms,tmp[[i]])}}
     taper_formula <- append_rhs.formula(~.,taper.terms)
   }else{if(taper.terms=="all"){
     taper.terms <- list_rhs.formula(formula)
     taper_formula <- append_rhs.formula(~.,taper.terms)
   }else{
    if(!inherits(taper.terms,"formula")){
      stop('taper.terms must be "dependent", "all" or a formula of terms.')
    }
    taper.terms <- list_rhs.formula(taper.terms)
    taper_formula <- append_rhs.formula(~.,taper.terms)
   }}
  }else{
    taper_formula <- taper.terms
    taper.terms <- list_rhs.formula(taper.terms)
  }
  taper.stats <- summary(append_rhs.formula(nw ~.,taper.terms))

# if(is.logical(taper.terms)){
#  if(length(taper.terms)!=length(ostats)){stop("The length of taper.terms must match that of the list of terms.")}
# }

  # set tapering coefficient
  tau <- switch(family,
    "stereo"={
      if(is.null(beta) & is.null(tau)){
        1
      }else{
       if(is.null(tau)){
        beta
       }else{
        tau
       }
      }},
      {if(is.null(tau)){
      if(is.null(beta)){
        1 / (r^2 * pmax(1,abs(taper.stats)))
      }else{
        1 / beta^2
      }}else{tau}}
  )
  tau <- tau*tau
  npar <- length(ostats)
  names(tau) <- paste0("M4(",names(taper.stats),")")
  
  taper_terms <- switch(family,
    "stereo"=statnet.common::nonsimp_update.formula(taper_formula,.~Stereo(~.,coef=.taper.coef,m=.taper.center),
              environment(), from.new=TRUE), 
             statnet.common::nonsimp_update.formula(taper_formula,.~M4(~.,coef=.taper.coef,m=.taper.center),
              environment(), from.new=TRUE) 
      )

  if(length(list_rhs.formula(formula))==length(taper.terms)){
    newformula <- switch(family,
      "stereo"=statnet.common::nonsimp_update.formula(formula,.~Stereo(~.,coef=.taper.coef,m=.taper.center),
                environment(), from.new=TRUE), 
               statnet.common::nonsimp_update.formula(formula,.~M4(~.,coef=.taper.coef,m=.taper.center),
                environment(), from.new=TRUE) 
	       )
#   newformula <- append_rhs.formula(newformula,taper_terms, environment()) 
  }else{
    trimmed_formula=suppressWarnings(filter_rhs.formula(formula, function(term,t.terms){all(term != taper.terms)}, taper.terms))
#   newformula <- append_rhs.formula(trimmed_formula,taper_terms, environment()) 
    newformula <- append_rhs.formula(formula,taper_terms, environment()) 
  }
  newformula <- append_rhs.formula(formula,taper_terms, environment()) 

  env <- new.env(parent=environment(formula))
  env$.taper.center <- taper.stats
  env$.taper.coef <- tau
  environment(newformula) <- env

  message(sprintf("The tapering formula is:\n %s", paste(deparse(newformula), sep="\n", collapse = "\n")))
  message("The (natural) tapering parameters are:")
  for(i in seq_along(tau)){
    message(sprintf(" %s : %f",names(tau)[i],tau[i]))
  }
  message("\n")
  
  control$MCMLE.metric <- "m4"
  if(control$MCMLE.termination == "Hotelling") control$MCMLE.termination <- "confidence"

  if(is.null(control$init)){
    fit0 <- c(ergm(formula, control=control, estimate="MPLE")$coef,tau)
    if(length(tau)==length(fit0)-npar-1){
     names(fit0)[(npar+1):length(fit0)] <- names(tau)
    }
    control$init <- fit0
  }

  # fit ergm
  if(is.null(target.stats)){
    fit <- ergm(newformula, control=control, ...)
  }else{
    fit <- ergm(newformula, control=control, target.stats=ostats, offset.coef=tau, ...)
  }
  
# # post processs fit to alter Hessian etc
# sample <- fit$sample[[1]][,1:npar,drop=FALSE]
# if(is.null(tapering.centers)){
#   hesstau <- ostats - ostats
##  hesstau[names(ostats) %in% names(tau)] <- tau
#   nm <- match(names(ostats),names(tau))
#   hesstau[seq_along(hesstau)[!is.na(nm)]] <- tau[nm[!is.na(nm)]]
#   hess <- .tapered.hessian(sample, hesstau)
#   if(is.curved(fit)){
#     curved_m <- ergm_model(newformula)
#     curved_m <- .tapered.curved.hessian(hess,fit$coef,curved_m$etamap)
#     fit$hessian[colnames(fit$hessian) %in% colnames(curved_m),rownames(fit$hessian) %in% rownames(curved_m)] <- curved_m
#   }else{
#     fit$hessian <- hess
#   }
    fit$covar <- -MASS::ginv(fit$hessian)
# }
  fit$tapering.centers <- ostats[as.character(unlist(taper.terms))]
  fit$tapering.coef <- tau
  fit$orig.formula <- formula
  class(fit) <- c("tapered.ergm",family,class(fit))
  
  fit
}
