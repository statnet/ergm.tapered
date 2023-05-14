#' Simulates from a Tapered ERGM formula
#' @param formula An ergm formula to fit
#' @param coef Vector of parameter values for the model from which the
#'   sample is to be drawn.
#' @param nsim Number of networks to be randomly drawn from the given
#' distribution on the set of all networks, returned by the Metropolis-Hastings
#' algorithm.
#' @param r The scaling factor to use for the heuristic of setting beta equal to r standard deviations of the observed statistics
#' @param beta The tapering parameters, expressed as in Fellows and Handcock (2017). If not NULL, these override the heuristics (r).
#' @param tau The tapering parameters, expressed as natural parameters. If not NULL, these override the beta and the heuristics (r).
#' @param tapering.centers The centers of the tapering terms. If null, these are taken to be the mean value parameters.
#' Equivalently, this vector is the mean-value parameter values for the
#' model.  If this is given, the algorithm finds the natural
#' parameter values corresponding to these mean-value parameters.
#' If \code{NULL}, the mean-value parameters used are the observed
#' statistics of the network in the formula.
#' }
#' @param family The type of tapering used. This should either be the \code{stereo} or \code{taper}, the 
#' tapering model of Fellows and Handcock (2017).
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
#' @param seed the random seed to use; see [simulate()].
#' @param control An object of class control.ergm. Passed to the ergm function.
#' @param verbose A `logical`: if this is
#' \code{TRUE}, the program will print out additional
#' information about the progress of estimation.
#' @param ... Additional arguments to \code{\link{ergm}}.
#' @returns
#' An object of class c('tapered.ergm','ergm') containing the fit model. In addition to all of the ergm items, 
#' this object contains tapering.centers, tapering.coef and orig.formula. tapering.centers are the centers for the tapering term.
#' tapering.coef are the tapering coefficients = 1/ beta^2. orig.formula is the formula passed into ergm.tapered.
#' @importFrom stats var as.formula
#' @references \itemize{ 
#' * Fellows, I. and M. S. Handcock (2017), 
#' Removing Phase Transitions from Gibbs Measures. Volume 54 of 
#' Proceedings of Machine Learning Research, Fort Lauderdale,
#' FL, USA, pp. 289–297. PMLR.
#' * Blackburn, B. and M. S. Handcock (2022), 
#' Practical Network Modeling via Tapered Exponential-family Random Graph Models.
#' Journal of Computational and Graphical Statistics
#' \doi{10.1080/10618600.2022.2116444}.
#' 
#' }
#' @examples 
#' \dontrun{
#' data(sampson)
#' fit <- ergm.tapered(samplike ~ edges + triangles())
#' summary(fit)
#' }
#' @export
simulate_ergm.tapered <- function(formula, coef, nsim=1, r=2, beta=NULL, tau=NULL, tapering.centers=NULL,
			 family="taper", taper.terms="all",
                         response=NULL, constraints=~., reference=~Bernoulli, seed=NULL,
                         control = control.simulate.formula(), verbose=FALSE, ...){

  # Needed by ergm.getMCMCsample
  match.llik.arg.pars <- c("MCMLE.metric","MCMLE.termination","MCMC.esteq.exclude.statistics",STATIC_TAPERING_CONTROLS)
  for(arg in match.llik.arg.pars)
    control$loglik[arg]<-list(control[[arg]])

  # Determine the dyadic independence terms
  nw <- ergm.getnetwork(formula)
  m<-ergm_model(formula, nw, response=response, ...)

  # set tapering terms
  if(is.character(taper.terms) & length(taper.terms)==1){
   if(taper.terms=="dependent"){
     a <- sapply(m$terms, function(term){is.null(term$dependence) || term$dependence})
     taper.terms <- list_rhs.formula(formula)
     tmp <- taper.terms
     taper.terms <- NULL
     for(i in seq_along(tmp)){if(a[i]){taper.terms <- c(taper.terms,tmp[[i]])}}
     taper_formula <- append_rhs.formula(~.,taper.terms, env=NULL)
     trimmed_formula=suppressWarnings(filter_rhs.formula(formula, function(term,taper.terms){all(term != taper.terms)}, taper.terms))
     reformula <- append_rhs.formula(trimmed_formula,taper_formula, env=NULL) 
   }else{if(taper.terms=="all"){
     taper.terms <- list_rhs.formula(formula)
     taper_formula <- append_rhs.formula(~.,taper.terms, env=NULL)
     trimmed_formula=suppressWarnings(filter_rhs.formula(formula, function(term,taper.terms){all(term != taper.terms)}, taper.terms))
     reformula <- formula
   }else{
    if(!inherits(taper.terms,"formula")){
      stop('taper.terms must be "dependent", "all" or a formula of terms.')
    }
    taper.terms <- list_rhs.formula(taper.terms)
    taper_formula <- append_rhs.formula(~.,taper.terms, env=NULL)
    trimmed_formula=suppressWarnings(filter_rhs.formula(formula, function(term,taper.terms){all(term != taper.terms)}, taper.terms))
    reformula <- append_rhs.formula(trimmed_formula,taper_formula, env=NULL) 
   }}
  }else{
    taper_formula <- taper.terms
    taper.terms <- list_rhs.formula(taper.terms)
    trimmed_formula=suppressWarnings(filter_rhs.formula(formula, function(term,taper.terms){all(term != taper.terms)}, taper.terms))
    reformula <- append_rhs.formula(trimmed_formula,taper_formula, env=NULL) 
  }
  attr(taper.terms,"env") <- NULL
  if(is.null(tapering.centers)){
    taper.stats <- summary(append_rhs.formula(nw ~.,taper.terms), response=response, ...)
  }else{
    taper.stats <- tapering.centers
  }

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
  if(family=="stereo"){
    names(tau) <- "stereo"
  }else{
    names(tau) <- names(taper.stats)
  }
  
  fixed <- TRUE
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
    newformula <- append_rhs.formula(trimmed_formula,taper_terms, env=NULL) 
  }
  ostats <- summary(reformula, response=response, ...)
  if(!is.null(tapering.centers)){
    tmp <- match(names(ostats), names(tapering.centers))
    if(any(is.na(tmp))){
      stop(paste("tapering.centers needs to have a named center for each statistic in the model:",
           names(ostats)))
    }
    ostats <- tapering.centers[tmp]
  }
  npar <- length(ostats)

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
  
  re.names <- names(summary(newformula))

  # fit ergm
  fit <- simulate(newformula, nsim=nsim, seed=seed, coef=coef, control=control,
              response=response, constraints=constraints, reference=reference, verbose=verbose, ...)
  
  if(is.matrix(fit)){
    a <- grep("Taper_Penalty",colnames(fit),fixed=TRUE)
    if(length(a) > 0){
      fit0 <- fit[,-a,drop=FALSE]
      attr(fit0, "monitored") <- attr(fit0, "monitored")[-a]
      attr(fit0, "formula") <- attr(fit0, "formula")
      attr(fit0, "constraints") <- attr(fit0, "constraints")
      attr(fit0, "reference") <- attr(fit0, "reference")
      fit <- fit0
    }
  }
  attr(fit, "tapering.centers") <- taper.stats
  attr(fit, "tapering.centers.o") <- ostats
  attr(fit, "tapering.centers.t") <- taper.terms
  attr(fit, "tapering.coef") <- tau
  attr(fit, "r") <- r

  fit
}
