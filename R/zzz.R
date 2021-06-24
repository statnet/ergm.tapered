.onUnload <- function(libpath){
  library.dynam.unload("ergm.tapered",libpath)
}
.onLoad <- function(libname, pkgname){
  utils::globalVariables(names=c(".","NO_LOGLIK_MESSAGE","NO_NULL_IMPLICATION"))
}
