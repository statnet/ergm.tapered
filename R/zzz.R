.onUnload <- function(libpath){
  library.dynam.unload("ergm.tapered",libpath)
}
.onAttach <- function(lib, pkg){
  if("package:ergm" %in% search()){
    packageStartupMessage("Detaching prior version of 'ergm'.")
    detach("package:ergm", unload=TRUE)
  }
}
