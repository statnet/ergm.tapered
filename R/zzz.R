.onUnload <- function(libpath){
  library.dynam.unload("ergm.tapered",libpath)
}

.onLoad <- function(libname, pkgname){
  utils::globalVariables(names=c(".","NO_LOGLIK_MESSAGE","NO_NULL_IMPLICATION"))
  .RegisterKeywords()
  if(utils::packageVersion("ergm") != "4.3.6983"){
    stop("ergm.tapered requires a variant of the current version of 'ergm' (the 'tapered' branch).\n  To install it use:\n devtools::install_github('statnet/ergm', ref='tapered')")
  }
}


.RegisterKeywords <- function() {
  ergm_keyword(name="tapered", short="taper", description="applies tapering to specified terms", popular=FALSE, package="ergm.tapered")
}
