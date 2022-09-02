.onUnload <- function(libpath){
  library.dynam.unload("ergm.tapered",libpath)
}


.onLoad <- function(libname, pkgname){
  utils::globalVariables(names=c(".","NO_LOGLIK_MESSAGE","NO_NULL_IMPLICATION"))
  .RegisterKeywords()
}


.RegisterKeywords <- function() {
  ergm_keyword(name="tapered", short="taper", description="applies tapering to specified terms", popular=FALSE, package="ergm.tapered")
}
