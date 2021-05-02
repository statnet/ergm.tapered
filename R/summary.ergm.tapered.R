#' Summarizing Tapered and Kurtosis ERGM Model Fits
#' 
#' \code{\link[base]{summary}} method for [`ergm.tapered`] fits.
#' 
#' \code{\link{summary.ergm.tapered}} tries to be smart about formatting the
#' coefficients, standard errors, etc.
#' 
#' @aliases print.summary.ergm.tapered
#' @param object an object of class \code{"ergm.tapered"}, usually, a result
#'   of a call to \code{\link{ergm.tapered}} or \code{\link{ergm.kurtosis}}.
#' @param digits Significant digits for coefficients
#' @param correlation logical; if \code{TRUE}, the correlation matrix
#'   of the estimated parameters is returned and printed.
#' @param covariance logical; if \code{TRUE}, the covariance matrix of
#'   the estimated parameters is returned and printed.
#' @param total.variation logical; if \code{TRUE}, the standard errors
#'   reported in the \code{Std. Error} column are based on the sum of
#'   the likelihood variation and the MCMC variation. If \code{FALSE}
#'   only the likelihood varuation is used. The \eqn{p}-values are
#'   based on this source of variation.
#' @param \dots Arguments to \code{\link{logLik.ergm}}
#' @return The function \code{\link{summary.ergm}} computes and
#'   returns a list of summary statistics of the fitted
#'   \code{\link{ergm}} model given in \code{object}. Note that for
#'   backwards compatibility, it returns two coefficient tables:
#'   `$coefs` which does not contain the z-statistics and
#'   `$coefficeints` which does (and is therefore more similar to
#'   those returned by [summary.lm()]).
#' @seealso network, ergm, print.ergm.  The model fitting function
#'   \code{\link{ergm}}, \code{\link{summary}}.
#' 
#' Function \code{\link{coef}} will extract the matrix of coefficients with
#' standard errors, t-statistics and p-values.
#' @keywords regression models
#' @examples
#' 
#'  data(florentine)
#' 
#'  x <- ergm(flomarriage ~ density)
#'  summary(x)
#' 
#' @export
summary.ergm.tapered <- function (object, ..., 
                          correlation=FALSE, covariance=FALSE,
                          total.variation=TRUE)
{
  if("digits" %in% names(list(...))) warn("summary.ergm() no longer takes a digits= argument.")
  
  control <- object$control
  pseudolikelihood <- object$estimate=="MPLE"
  independence <- NVL(object$MPLE_is_MLE, is.dyad.independent(object))
  
  if(any(is.na(object$coef)) & !is.null(object$mplefit)){
     object$coef[is.na(object$coef)] <-
     object$mplefit$coef[is.na(object$coef)]
  }

  nodes<- network.size(object$network)

  cnames.all <- param_names(object)

  cnames <- grep("Taper_Penalty",cnames.all,fixed=TRUE,invert=TRUE)

  ans <- list(formula=object$formula,
              correlation=correlation,
              degeneracy.value = object$degeneracy.value,
              offset = object$offset[cnames],
              drop = NVL(object$drop[cnames], rep(0,length(object$offset[cnames]))),
              estimable = NVL(object$estimable[cnames], rep(TRUE,length(object$offset[cnames]))),
              covariance=covariance,
              pseudolikelihood=pseudolikelihood,
              independence=independence,
              estimate=object$estimate,
              control=object$control)
  
  ans$samplesize <- switch(object$estimate,
                           EGMME = NVL3(control$EGMME.main.method, switch(.,
                             `Gradient-Descent`=control$SA.phase3n,
                             stop("Unknown estimation method. This is a bug."))),
                           MPLE = NA,
                           CD=,
                           MLE = NVL3(control$main.method, switch(.,
                             CD=control$MCMC.samplesize,
                             `Stochastic-Approximation`=,
                               MCMLE=control$MCMC.samplesize,
                             `Robbins-Monro`=control$RM.phase3n,
                             `Stepping`=control$Step.MCMC.samplesize,
                             stop("Unknown estimation method. This is a bug."))),
                           stop("Unknown estimate type. This is a bug.")
                           )
                              

  ans$iterations <- switch(object$estimate,
                           EGMME = NVL3(control$EGMME.main.method, switch(.,
                             `Gradient-Descent`=NA,
                             stop("Unknown estimation method. This is a bug."))),
                           MPLE = NA,
                           CD=control$CD.maxit,
                           MLE = NVL3(control$main.method, switch(.,
                               `Stochastic-Approximation`=NA,
                             MCMLE=paste(object$iterations, "out of", control$MCMLE.maxit),
                             CD=control$CD.maxit,
                             `Robbins-Monro`=NA,
                             `Stepping`=NA,
                             stop("Unknown estimation method. This is a bug."))),
                           stop("Unknown estimate type. This is a bug.")
                           )
  
  nodes<- network.size(object$network)
  dyads<- sum(as.rlebdm(object$constrained, object$constrained.obs, which="informative"))
  df <- length(object$coef[cnames])

  asycov <- vcov(object, sources=if(total.variation) "all" else "model")[cnames,cnames]
# asycov <- vcov(object, sources="estimation")[cnames,cnames]
# asycov <- object$covar
# asycov <- MASS::ginv(-object$covar[cnames,cnames])
  asyse <- sqrt(diag(asycov))
  # Convert to % error  
  est.se <- sqrt(diag(vcov(object, sources="estimation")))[cnames]
  mod.se <- sqrt(diag(vcov(object, sources="model")))[cnames]
# mod.se <- asyse
  tot.se <- sqrt(diag(vcov(object, sources="all")))[cnames]
# tot.se <- sqrt(est.se^2 + mod.se^2)
  est.pct <- rep(NA,length(est.se))
  if(any(!is.na(est.se))){
    # We want (sqrt(V.model + V.MCMC)-sqrt(V.model))/sqrt(V.model + V.MCMC) * 100%,
    est.pct[!is.na(est.se)] <- ifelse(est.se[!is.na(est.se)]>0, round(100*(tot.se[!is.na(est.se)]-mod.se[!is.na(est.se)])/tot.se[!is.na(est.se)]), 0)
  }

  rdf <- dyads - df
  zval <- object$coef[cnames] / asyse
  pval <- 2 * pnorm(q=abs(zval), lower.tail=FALSE)
  
  count <- 1
  coefmat <- cbind(
    `Estimate` = coef(object)[cnames],
    `Std. Error` = asyse,
    `MCMC %` = est.pct,
    `z value` = zval,
    `Pr(>|z|)` = pval)

  devtext <- "Deviance:"
  if (object$estimate!="MPLE" || !independence || object$reference != as.formula(~Bernoulli)) {
    if (pseudolikelihood) {
      devtext <- "Pseudo-deviance:"
      ans$message <- "\nWarning:  The standard errors are based on naive pseudolikelihood and are suspect.\n"
    } 
    else if(object$estimate == "MLE" && any(is.na(est.se) & !ans$offset & !ans$drop==0 & !ans$estimable) && 
                      (!independence || control$force.main) ) {
      ans$message <- "\nWarning:  The standard errors are suspect due to possible poor convergence.\n"
    }
  } else {
    ans$message <- "\nFor this model, the pseudolikelihood is the same as the likelihood.\n"
  }
  null.lik<-try(logLikNull(object,...), silent=TRUE)
  if(inherits(null.lik,"try-error")){
    null.lik<-object$null.lik
  }
  mle.lik<-try(logLik(object,...), silent=TRUE)
  if(inherits(mle.lik,"try-error")){
    mle.lik<-null.lik+object$loglikelihood
  }

  ans$null.lik.0 <- is.na(null.lik)

  if(!inherits(mle.lik,"try-error")){

    ans$devtable <- matrix(c(if(is.na(null.lik)) 0 else -2*null.lik, -2*mle.lik,
                             c(dyads, rdf)), 2,2, dimnames=list(c("Null","Residual"),
                                                                c("Resid. Dev", "Resid. Df")))
    ans$devtext <- devtext
        
    ans$aic <- AIC(mle.lik)
    ans$bic <- BIC(mle.lik)
    ans$mle.lik <- ERRVL(mle.lik, NA)
    ans$null.lik <- ERRVL(null.lik, NA)
  }else ans$objname<-deparse(substitute(object))

  ans$coefs <- as.data.frame(coefmat)[cnames,-3] # For backwards compatibility.
  ans$coefficients <- as.data.frame(coefmat)
  ans$asycov <- asycov
  ans$asyse <- asyse
  ans$orig.formula <- object$orig.formula
  ans$r <- object$r
  a <- grep("Taper_Penalty",cnames.all,fixed=TRUE)
  if(length(a)>0){
    ans$Taper_Penalty <- object$coef[a]
    ans$r <- object$r/sqrt(3*exp(log(2)*(ans$Taper_Penalty-2))/(1+exp(log(2)*(ans$Taper_Penalty-2))))
  }
  a <- grep("Var(",cnames.all,fixed=TRUE)
  if(length(a)>0){
    a <- -1/(object$coef[a]*pmax(1,object$tapering.centers))
    a[is.nan(a) | a < 0] <- 0
    ans$r <- mean(sqrt(a))
  }
  class(ans) <- "summary.ergm.tapered"
  ans
}
