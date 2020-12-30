InitWtErgmTerm.Stereo <- function(nw, arglist, response=NULL, ...){
  out <- InitErgmTerm.Stereo(nw, arglist, response=response, ...)
  out$name <- "wttaper_term"
  out
}
