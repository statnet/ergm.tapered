#' Fits a Mean-Valued Tapered ERGM
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
ergm.mvtapered <- function(formula, r=2, mv=NULL, tapering.centers=NULL, 
			 family="taper", taper.terms="dependent",
                         control = control.ergm(), ...){
  
  # Determine the dyadic independence terms
  nw <- ergm.getnetwork(formula)
  m<-ergm_model(formula, nw)
  if(is.character(taper.terms)){
   taper.terms <- switch(taper.terms,
    "dependent"={
     sapply(m$terms, function(term) is.null(term$dependence) || term$dependence)
    },
     rep(TRUE,length(m$terms))
   )
  }
# if(length(tapering.centers) != length(taper.terms)){
#  stop("'tapering.centers' is the wrong length.")
# }

  control$SAN.maxit <- 0
  control$MCMLE.termination <- "Hotelling"

  if(is.null(tapering.centers))
    ostats <- summary(formula)
  else
    ostats <- tapering.centers
  
  # set tapering coefficient
  target.stats <- switch(family,
    "stereo"={
      if(is.null(beta)){
        1
      }else{
        beta
      }},
      {if(length(mv)==1 & is.numeric(mv)){
        c(ostats,mv*abs(ostats[taper.terms]))
      }else{
        c(ostats, mv)
      }}
  )
  npar <- length(ostats)
  coef <- rep(1,length(ostats))
  
  terms <- list_rhs.formula(formula)
  terms[!taper.terms] <- NULL
  attr(terms,"sign") <- attr(terms,"sign")[taper.terms]
  taper_formula=append_rhs.formula(~.,terms)

  # do some formula magic

  taper_formula <- switch(family,
    "stereo"=statnet.common::nonsimp_update.formula(formula,
              .~Stereo(~.,coef=.taper.coef,m=.taper.center), environment(), 
              from.new=TRUE),
             statnet.common::nonsimp_update.formula(taper_formula,
              ~ Var(~.,coef=.taper.coef,m=.taper.center), environment(), 
              from.new=TRUE) 
  )
  newformula <- append_rhs.formula(formula, taper_formula)
  env <- new.env(parent=environment(formula))
  env$.taper.center <- ostats[taper.terms]
  env$.taper.coef <- coef[taper.terms]
  environment(newformula) <- env
  
# target.stats=c(ostats,1*ostats)
  names(target.stats) <- names(summary(newformula))
  
  # fit ergm
  fit <- ergm(newformula, target.stats=target.stats, control=control,...)
  
  names(coef) <- names(ostats)
# # post processs fit to alter Hessian etc
# sample <- fit$sample[[1]][,1:npar,drop=FALSE]
# if(is.null(tapering.centers)){
#   hess <- .taperedMV.hessian(sample, coef)
#   if(is.curved(fit)){
#     curved_m <- ergm_model(formula)
#     curved_m <- .taperedMV.curved.hessian(hess,fit$coef,curved_m$etamap)
#     fit$hessian[colnames(fit$hessian) %in% colnames(curved_m),rownames(fit$hessian) %in% rownames(curved_m)] <- curved_m
#   }else{
#     fit$hessian <- hess
#   }
#   fit$covar <- -MASS::ginv(fit$hessian)
# }
  fit$tapering.centers <- ostats[taper.terms]
  fit$tapering.coef <- coef[taper.terms]
  fit$orig.formula <- formula
  class(fit) <- c("tapered.ergm",family,class(fit))
  
  fit
}

.taperedMV.hessian <- function(sample, coef){
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

.taperedMV.curved.hessian <- function(hess,theta,etamap){
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
