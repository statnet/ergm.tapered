#' Fits a Tapered ERGM
#' @param formula An ergm formula to fit
#' @param r The scaling factor to use for the hueristic of setting beta equal to r standard deviations of the observed statistics
#' @param beta The tapering parameters, expressed as in Fellows and Handcock (2017). If not NULL, these override the hueristics (r).
#' @param tau The tapering parameters, expressed as natural parameters. If not NULL, these override the beta and the hueristics (r).
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
#' @param family The type of tapering used. This should either be the \code{stereo} or \code{taper}, the 
#' tapering model of Fellows and Handcock (2016).
#' @param taper.terms Specification of the tapering used. If the character variable "dependence" then all the dependent
#' terms are tapered. If the character variable "all" then all terms are tapered.
#' It can also be the RHS of a formula giving the terms to be tapered. 
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
ergm.tapered <- function(formula, r=2, beta=NULL, tau=NULL, tapering.centers=NULL, target.stats=NULL,
			 family="taper", taper.terms="dependent",
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
     for(i in seq_along(taper.terms)){if(!a[i]){taper.terms[[i]] <- NULL}}
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
  npar <- length(ostats)
  names(tau) <- names(taper.stats)
  
  taper_terms <- switch(family,
    "stereo"=statnet.common::nonsimp_update.formula(taper_formula,.~Stereo(~.,coef=.taper.coef,m=.taper.center),
              environment(), from.new=TRUE), 
             statnet.common::nonsimp_update.formula(taper_formula,.~Taper(~.,coef=.taper.coef,m=.taper.center),
              environment(), from.new=TRUE) 
      )
# trimmed_formula=filter_rhs.formula(formula, function(term,t.terms){print(length(term));print(term);out <- FALSE;for(i in seq_along(t.terms)){out <- out | (paste(term,collapse="") == paste(taper.terms[[i]],collapse=""))};print(out);out},
#                                    taper.terms)

  if(length(list_rhs.formula(formula))==length(taper.terms)){
    newformula <- switch(family,
      "stereo"=statnet.common::nonsimp_update.formula(formula,.~Stereo(~.,coef=.taper.coef,m=.taper.center),
                environment(), from.new=TRUE), 
               statnet.common::nonsimp_update.formula(formula,.~Taper(~.,coef=.taper.coef,m=.taper.center),
                environment(), from.new=TRUE) 
	       )

  }else{
    trimmed_formula=suppressWarnings(filter_rhs.formula(formula, function(term,t.terms){all(term != taper.terms)}, taper.terms))
    newformula <- append_rhs.formula(trimmed_formula,taper_terms, environment()) 
  }

  env <- new.env(parent=environment(formula))
  env$.taper.center <- taper.stats
  env$.taper.coef <- tau
  environment(newformula) <- env

  message("The (natural) tapering parameters are:")
  for(i in seq_along(tau)){
    message(sprintf(" %s : %f",names(tau)[i],tau[i]))
  }
  
  if(control$MCMLE.termination == "Hotelling") control$MCMLE.termination <- "confidence"

  print(newformula)
  # fit ergm
  if(is.null(target.stats)){
    fit <- ergm(newformula, control=control, ...)
  }else{
    fit <- ergm(newformula, control=control, target.stats=ostats, offset.coef=tau, ...)
  }
  
  
  # post processs fit to alter Hessian etc
  sample <- fit$sample[[1]][,1:npar,drop=FALSE]
  if(is.null(tapering.centers)){
    hess <- .tapered.hessian(sample, tau)
    if(is.curved(fit)){
      curved_m <- ergm_model(formula)
      curved_m <- .tapered.curved.hessian(hess,fit$coef,curved_m$etamap)
      fit$hessian[colnames(fit$hessian) %in% colnames(curved_m),rownames(fit$hessian) %in% rownames(curved_m)] <- curved_m
    }else{
      fit$hessian <- hess
    }
    fit$covar <- -MASS::ginv(fit$hessian)
  }
  fit$tapering.centers <- ostats[taper.terms]
  fit$tapering.coef <- tau
  fit$orig.formula <- formula
  class(fit) <- c("tapered.ergm",family,class(fit))
  
  fit
}

.tapered.hessian <- function(sample, coef){
  cv <- var(sample)
  
  np <-ncol(cv)
  B <- sweep(cv, 2, 2*coef, "*")
  I <- diag(rep(1,np))
  inv <- MASS::ginv(B+I)
  
  #derivative of mean value parameters
  dmu <- inv %*% cv
  
  #second derivative of log likelihoods
  ddll <- diag(rep(0,np))
  for(i in 1:np){
    for(j in 1:np){
      ddll[i,j] <- -dmu[i,j] - sum(2*coef*dmu[,i]*dmu[,j]) 
    }
  }
  ddll
}

.tapered.curved.hessian <- function(hess,theta,etamap){
  eta <- rep(0,etamap$etalength)
  ec <- etamap$canonical
  eta[ec[ec>0]] <- theta[ec>0]
  nstats <- sum(ec>0) + length(etamap$curved)

  if(length(etamap$curved)>0) {
    ch <- matrix(0, nrow=length(eta), ncol=nstats)
    namesch <- rep("", nstats)
    icurved <- 0
    istat <- 0
    for(i in 1:nstats) {
      istat <- istat + 1
      if(ec[istat] > 0){
       ch[istat,i] <- 1
       namesch[i] <- names(theta)[istat]
      }else{
       icurved <- icurved + 1
       cm <- etamap$curved[[icurved]]
       ch[cm$to,i] <- cm$map(theta[cm$from],length(cm$to),cm$cov)
#      (un)scale by linear coefficient
       ch[cm$to,i] <- ch[cm$to,i] / theta[cm$from][1]
       namesch[i] <- names(theta)[istat]
       istat <- istat + length(cm$from) - 1
      }
    }
    ch <- t(ch) %*% hess %*% ch
    colnames(ch) <- namesch
    rownames(ch) <- namesch
  }else{
    ch <- hess
  }
  ch
}
