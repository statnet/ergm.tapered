.onUnload <- function(libpath){
  library.dynam.unload("ergm.tapered",libpath)
}
.onAttach <- function(lib, pkg){
  if("package:ergm" %in% search()){
    message("Detaching prior version of 'ergm'.")
    detach("package:ergm", unload=TRUE)
  }
}
