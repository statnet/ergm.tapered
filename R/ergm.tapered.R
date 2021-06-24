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
#' @param response {Name of the edge attribute whose value is to be
#' modeled in the valued ERGM framework. Defaults to \code{NULL} for
#' simple presence or absence, modeled via a binary ERGM.}
#' @param reference {A one-sided formula specifying
#' the reference measure (\eqn{h(y)}) to be used.
#' See help for \link[=ergm-references]{ERGM reference measures} implemented in the
#' \strong{\link[=ergm-package]{ergm}} package.}
#' @param constraints {A formula specifying one or more constraints
#' on the support of the distribution of the networks being modeled,
#' using syntax similar to the \code{formula} argument, on the
#' right-hand side. See \link[=ergm-constraints]{ERGM constraints} 
#' \code{\link{ergm}} for details.}
#' @param control An object of class control.ergm. Passed to the ergm function.
#' @param fixed A `logical`: if this is \code{TRUE}, the tapering is fixed at the passed values. 
#' If this is \code{FALSE}, the tapering is estimated using a kurtosis penalized likelihood.
#' @param verbose A `logical`: if this is
#' \code{TRUE}, the program will print out additional
#' information about the progress of estimation.
#' @param eval.loglik {Logical:  If TRUE, evaluate the log-likelihood associated with the fit.}
#' @param ... Additional arguments to \code{\link{ergm}}.
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
			 family="taper", taper.terms="all",
                         response=NULL, constraints=~., reference=~Bernoulli,
                         control = control.ergm.tapered(), fixed=TRUE, verbose=FALSE, eval.loglik=TRUE, ...){

  if(!methods::is(control,"control.ergm.tapered")){
    stop("The control variable must be of class 'control.ergm.tapered'.")
  }
  # Enforce Kpenalty metric
  if(!fixed){
    otaper.terms <- taper.terms
    control$MCMLE.metric <- "Kpenalty"
  }
  control["MCMC.esteq.exclude.statistics"] <- "Taper_Penalty"
  # Statistics to exclude from the estimating equations
  # Needed by ergm.getMCMCsample
  match.llik.arg.pars <- c("MCMLE.metric","MCMC.esteq.exclude.statistics",STATIC_TAPERING_CONTROLS)
  for(arg in match.llik.arg.pars)
    control$loglik[arg]<-list(control[[arg]])

  # Determine the dyadic independence terms
  nw <- ergm.getnetwork(formula)
  m<-ergm_model(formula, nw, response=response, ...)

  if(is.null(tapering.centers)) tapering.centers <- target.stats

  if(is.null(tapering.centers))
    ostats <- summary(formula, response=response, ...)
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
  taper.stats <- summary(append_rhs.formula(nw ~.,taper.terms), response=response, ...)

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
  if(family=="stereo"){
    names(tau) <- "stereo"
  }else{
    names(tau) <- names(taper.stats)
  }
  
  if(!fixed & is.null(control$init)){
    tcontrol <- control
    tcontrol$MCMLE.metric <- "lognormal"
    tcontrol$MCMLE.maxit <- 10
    tcontrol$MCMLE.steplength <- 0.5
    tcontrol.orig <- tcontrol$MCMC.effectiveSize
#   tcontrol$MCMC.effectiveSize <- 100
    tcontrol$MCMC.effectiveSize.maxruns=8
    tfit <- ergm.tapered(formula, r=r, beta=beta, tau=tau, 
         family=family, taper.terms=otaper.terms,
         response=response, constraints=constraints, reference=reference,
         control = tcontrol, eval.loglik=FALSE, verbose=FALSE, ...)
    control$init <- tfit$coef
  }

  taper_terms <- switch(paste0(family,ifelse(fixed,"_fixed","_notfixed")),
    "stereo_fixed"=statnet.common::nonsimp_update.formula(taper_formula,.~Stereo(~.,coef=.taper.coef,m=.taper.center),
              environment(), from.new=TRUE), 
    "stereo_notfixed"=statnet.common::nonsimp_update.formula(taper_formula,.~Stereo(~.,coef=.taper.coef,m=.taper.center),
              environment(), from.new=TRUE), 
    "taper_fixed"=statnet.common::nonsimp_update.formula(taper_formula,.~Taper(~.,coef=.taper.coef,m=.taper.center),
              environment(), from.new=TRUE), 
    "taper_notfixed"=statnet.common::nonsimp_update.formula(taper_formula,.~Kpenalty(~.,coef=.taper.coef,m=.taper.center),
              environment(), from.new=TRUE),
             statnet.common::nonsimp_update.formula(taper_formula,.~Taper(~.,coef=.taper.coef,m=.taper.center),
              environment(), from.new=TRUE)
      )

  if(length(list_rhs.formula(formula))==length(taper.terms)){
    newformula <- switch(paste0(family,ifelse(fixed,"_fixed","_notfixed")),
      "stereo_fixed"=statnet.common::nonsimp_update.formula(formula,.~Stereo(~.,coef=.taper.coef,m=.taper.center),
                environment(), from.new=TRUE), 
      "stereo_notfixed"=statnet.common::nonsimp_update.formula(formula,.~Stereo(~.,coef=.taper.coef,m=.taper.center),
                environment(), from.new=TRUE), 
      "taper_fixed"=statnet.common::nonsimp_update.formula(formula,.~Taper(~.,coef=.taper.coef,m=.taper.center),
                environment(), from.new=TRUE), 
      "taper_notfixed"=statnet.common::nonsimp_update.formula(formula,.~Kpenalty(~.,coef=.taper.coef,m=.taper.center),
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

  if(verbose){
    message(sprintf("The tapering formula is:\n %s", paste(deparse(newformula), sep="\n", collapse = "\n")))
    message("The (natural) tapering parameters are:")
    for(i in seq_along(tau)){
      message(sprintf(" %s : %f",names(tau)[i],tau[i]))
    }
    message("\n")
  }
  
# control$main.hessian <- FALSE
# if(control$MCMLE.termination == "Hotelling") control$MCMLE.termination <- "confidence"
  re.names <- names(summary(newformula))
  if(!fixed){
    control$init <- c(control$init,1)
    names(control$init)[length(control$init)] <- "Taper_Penalty"
    control$init <- control$init[match(re.names,names(control$init))]
  }else{
    if(!is.null(control$init)){
      warning("check the names are ordered correctly (not coded yet)")
      print(control$init)
      re.names <- re.names[-length(re.names)]
      names(control$init) <- re.names
    }
  }

  # fit ergm
  if(is.null(target.stats)){
    fit <- ergm(newformula, control=control,
                response=response, constraints=constraints, reference=reference, eval.loglik=eval.loglik, verbose=verbose, ...)
  }else{
    fit <- ergm(newformula, control=control, target.stats=ostats, offset.coef=tau,
                response=response, constraints=constraints, reference=reference, eval.loglik=eval.loglik, verbose=verbose, ...)
  }
  
    # post processs fit to alter Hessian etc
    if(fixed){
      sample <- as.matrix(fit$sample)[,1:npar,drop=FALSE]
      fit$hessian <- fit$hessian[1:npar,1:npar]
      fcoef <- fit$coef[1:npar]
    }else{
      sample <- as.matrix(fit$sample)[,1:npar,drop=FALSE]
    # fit$hessian <- fit$hessian[1:npar,1:npar]
    # fit$hessian <- -fit$hessian
      fcoef <- fit$coef[1:npar]
    }
    colnames(sample) <- names(ostats)
    
    if(is.null(tapering.centers)){
      hesstau <- ostats - ostats
  #   hesstau[names(ostats) %in% names(tau)] <- tau
      nm <- match(names(ostats),names(tau))
      hesstau[seq_along(hesstau)[!is.na(nm)]] <- tau[nm[!is.na(nm)]]
      hess <- .tapered.hessian(sample, hesstau)
      if(is.curved(fit)){
        curved_m <- ergm_model(formula, nw, response=response, ...)
        curved_m <- .tapered.curved.hessian(hess,fcoef,curved_m$etamap)
        fit$hessian[colnames(fit$hessian) %in% colnames(curved_m),rownames(fit$hessian) %in% rownames(curved_m)] <- curved_m
        fit$covar[colnames(fit$hessian) %in% colnames(curved_m),rownames(fit$hessian) %in% rownames(curved_m)] <- MASS::ginv(curved_m)
      }else{
        fit$hessian <- hess
        fit$covar <- -MASS::ginv(fit$hessian)
      }
    }
 #}else{
 #    fit$covar <- -MASS::ginv(fit$hessian)
 #}
# fit$tapering.centers <- ostats[as.character(unlist(taper.terms))]

  fit$tapering.centers <- taper.stats
  fit$tapering.centers.o <- ostats
  fit$tapering.centers.t <- taper.terms
  fit$tapering.coef <- tau
  fit$r <- r
  fit$orig.formula <- formula
  class(fit) <- c("ergm.tapered",family,class(fit))
  
# fit$covar <- fit$est.cov
# fit$covar <- -MASS::ginv(fit$hessian)
# fit$est.cov <- fit$covar

  fit
}

.tapered.hessian <- function(sample, coef){
  cv <- var(sample)
  
  np <-ncol(cv)
  B <- sweep(cv, 2, 2*coef, "*")
  I <- diag(rep(1,np))
  inv <- MASS::ginv(B-I)
  
  #derivative of mean value parameters
  dmu <- inv %*% cv
  
  #second derivative of log likelihoods
  ddll <- diag(rep(0,np))
  colnames(ddll) <- colnames(sample)
  rownames(ddll) <- colnames(sample)
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
    exterms <- c(grep("Var(",namesch,fixed=TRUE), grep("Taper_Penalty",namesch,fixed=TRUE))
    if(length(exterms)>0){
     ch <- ch[-exterms,-exterms]
     namesch <- namesch[-exterms]
    }
    ch <- t(ch) %*% hess %*% ch
    colnames(ch) <- namesch
    rownames(ch) <- namesch
  }else{
    ch <- hess
  }
  ch
}
